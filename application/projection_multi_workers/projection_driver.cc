/*
 * Refactor the functions of original PhysBAM code "PROJECTION_DRIVER.cpp".
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include "nimbus_pcg_sparse_mpi.h" 
#include "projection_driver.h"
#include "projection_example.h"

using namespace PhysBAM;

namespace {
template<class TV> void Write_Substep_Helper(void* writer,
		const std::string& title, int substep, int level) {
	((PROJECTION_DRIVER<TV>*)writer)->Write_Substep(title,
			substep, level);
}
};

// A driver is initialized with a pinter to example class.
template<class TV> PROJECTION_DRIVER<TV>::PROJECTION_DRIVER(
		PROJECTION_EXAMPLE<TV>& example) :
	example(example) {
	DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this, &Write_Substep_Helper<TV>);
}

template<class TV> PROJECTION_DRIVER<TV>::~PROJECTION_DRIVER() {
	DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}

template<class TV> void PROJECTION_DRIVER<TV>::PrepareForProjection() {
	Initialize();
	projection_data = new ProjectionData<TV>;
	projection_internal_data = new ProjectionInternalData<TV>;
	T dt=example.Time_At_Frame(current_frame+1)-time;
	example.Set_Boundary_Conditions(time+dt);
	example.projection.p*=dt; // rescale pressure for guess
	example.boundary_scalar.Apply_Boundary_Condition_Face(example.mac_grid,
			example.face_velocities, time+dt);
	example.mpi_grid->Average_Common_Face_Data(example.face_velocities);
	example.projection.Compute_Divergence(
			typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP(example.face_velocities),
			example.projection.elliptic_solver);

	LAPLACE_UNIFORM<GRID<TV> >* laplace = example.projection.elliptic_solver;
	laplace->Find_Solution_Regions(); // flood fill
	typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
	typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
	typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
	typedef typename GRID<TV>::VECTOR_INT TV_INT;

	//ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array(laplace->number_of_regions);
	projection_data->matrix_index_to_cell_index_array
			= new ARRAY<ARRAY<TV_INT> >(laplace->number_of_regions);
	ARRAY<ARRAY<TV_INT> > &matrix_index_to_cell_index_array =
			*(projection_data->matrix_index_to_cell_index_array);

	T_ARRAYS_INT cell_index_to_matrix_index(laplace->grid.Domain_Indices(1));
	ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-1, laplace->number_of_regions);

	//ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array(laplace->number_of_regions);
	projection_data->A_array = new ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >(laplace->number_of_regions);
	ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > &A_array = *(projection_data->A_array);

	// ARRAY<VECTOR_ND<T> > b_array(laplace->number_of_regions);
	projection_data->b_array = new ARRAY<VECTOR_ND<T> >(laplace->number_of_regions);
	ARRAY<VECTOR_ND<T> > &b_array = *(projection_data->b_array);

	for (CELL_ITERATOR iterator(laplace->grid, 1); iterator.Valid(); iterator.Next())
		filled_region_cell_count(laplace->filled_region_colors(iterator.Cell_Index()))++;
	for (int color=1; color<=laplace->number_of_regions; color++)
		if (laplace->filled_region_touches_dirichlet(color)
				||laplace->solve_neumann_regions) {
			matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));
		}
	filled_region_cell_count.Fill(0); // reusing this array in order to make the indirection arrays
	DOMAIN_ITERATOR_THREADED_ALPHA<LAPLACE_UNIFORM<GRID<TV> >,TV>
			threaded_iterator(laplace->grid.Domain_Indices(1),
					laplace->thread_queue, 1, 1, 2, 1);

	// ARRAY<int,TV_INT> domain_index(laplace->grid.Domain_Indices(1),false);
	projection_data->domain_index = new ARRAY<int, TV_INT>(laplace->grid.Domain_Indices(1), false);
	ARRAY<int,TV_INT> &domain_index = *(projection_data->domain_index);

	for (int i=1; i<=threaded_iterator.domains.m; i++) {
		RANGE<TV_INT> interior_domain(threaded_iterator.domains(i));
		interior_domain.max_corner-=TV_INT::All_Ones_Vector();
		interior_domain.min_corner+=TV_INT::All_Ones_Vector();
		for (CELL_ITERATOR iterator(laplace->grid, interior_domain); iterator.Valid(); iterator.Next())
			domain_index(iterator.Cell_Index())=i;
	}
	ARRAY<ARRAY<INTERVAL<int> > > interior_indices(laplace->number_of_regions);
	ARRAY<ARRAY<ARRAY<INTERVAL<int> > > >
			ghost_indices(laplace->number_of_regions);
	for (int color=1; color<=laplace->number_of_regions; color++) {
		interior_indices(color).Resize(threaded_iterator.number_of_domains);
		ghost_indices(color).Resize(threaded_iterator.number_of_domains);
		for (int i=1; i<=threaded_iterator.domains.m; i++)
			ghost_indices(color)(i).Resize(2*TV::dimension);
	}
	laplace->laplace_mpi->Find_Matrix_Indices(filled_region_cell_count,
			cell_index_to_matrix_index, matrix_index_to_cell_index_array);
	RANGE<TV_INT> domain=laplace->grid.Domain_Indices(1);
	laplace->Find_A(domain, A_array, b_array, filled_region_cell_count,
			cell_index_to_matrix_index);
}

template<class TV> void PROJECTION_DRIVER<TV>::Initialize() {
	DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

	// setup time
	if (example.restart)
		current_frame=example.restart;
	else
		current_frame=example.first_frame;
	time=example.Time_At_Frame(current_frame);

	// mpi
	if (example.mpi_grid)
		example.mpi_grid->Initialize(example.domain_boundary);
	example.projection.elliptic_solver->mpi_grid=example.mpi_grid;
	if (example.mpi_grid)
		example.boundary=new BOUNDARY_MPI<GRID<TV>,T>(example.mpi_grid,example.boundary_scalar);
	else
		example.boundary=&example.boundary_scalar;

	//threading
	example.projection.elliptic_solver->thread_queue=example.thread_queue;

	// setup grids and velocities
	example.projection.Initialize_Grid(example.mac_grid);
	example.face_velocities.Resize(example.mac_grid);
	example.Initialize_Fields();

	// setup laplace
	example.projection.elliptic_solver->Set_Relative_Tolerance(1e-11);
	example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(200);
	example.projection.elliptic_solver->pcg.evolution_solver_type
			=krylov_solver_cg;
	example.projection.elliptic_solver->pcg.cg_restart_iterations=40;

	if (example.restart)
		example.Read_Output_Files(example.restart);

	// setup domain boundaries
	VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;
	constant_extrapolation.Fill(VECTOR<bool, 2>::Constant_Vector(true));
	example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
	example.Set_Boundary_Conditions(time); // get so CFL is correct

	if (!example.restart)
		Write_Output_Files(example.first_frame);
	output_number=example.first_frame;
}

/*
 // Original code to handle each region.
 template<class TV> void PROJECTION_DRIVER<TV>::
 LoopEachRegionSetUp() {
 LAPLACE_UNIFORM<GRID<TV> >* laplace = example.projection.elliptic_solver;
 for(int color=1;color<=laplace->number_of_regions;color++)
 if(filled_region_cell_count(color)>0
 && (laplace->filled_region_touches_dirichlet(color)||laplace->solve_neumann_regions)){
 printf("After enter the loop!\n");
 Solve_Subregion(
 laplace,
 interior_indices(color),
 ghost_indices(color),
 matrix_index_to_cell_index_array(color),
 A_array(color),
 b_array(color),
 color,
 &domain_index);
 }
 }
 */

template<class TV> void PROJECTION_DRIVER<TV>::ApplyPressureAndFinish() {
	typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
	LAPLACE_UNIFORM<GRID<TV> >* laplace = example.projection.elliptic_solver;
	T dt=example.Time_At_Frame(current_frame+1)-time;
	for (CELL_ITERATOR iterator(laplace->grid, 1); iterator.Valid(); iterator.Next()) {
		int filled_region_color=
				laplace->filled_region_colors(iterator.Cell_Index());
		if (filled_region_color>0
				&& !laplace->filled_region_touches_dirichlet(filled_region_color))
			laplace->u(iterator.Cell_Index())=0;
	}
	example.projection.Apply_Pressure(example.face_velocities, dt, time);
	example.projection.p/=dt; // rescale pressure for guess
	time+=dt;
	Write_Output_Files(++output_number);
	current_frame++;
}

template<class TV> void PROJECTION_DRIVER<TV>::PrepareForOneRegion() {
	int color = 1;
	LAPLACE_UNIFORM<GRID<TV> >* laplace = example.projection.elliptic_solver;
	laplace->pcg.Enforce_Compatibility( !laplace->filled_region_touches_dirichlet(color)
			&&laplace->enforce_compatibility);
	// Context.
	ARRAY<TV_INT> &matrix_index_to_cell_index =
			(*(projection_data->matrix_index_to_cell_index_array))(color);
	int number_of_unknowns=matrix_index_to_cell_index.m;
	SPARSE_MATRIX_FLAT_NXN<T> &A = (*(projection_data->A_array))(color);
	A.Negate();
	VECTOR_ND<T> &b = (*(projection_data->b_array))(color);
	b*=(T)-1;
	projection_data->x = new VECTOR_ND<T>(number_of_unknowns);
	VECTOR_ND<T> &x = *(projection_data->x);
	// Context done.
	VECTOR_ND<T> q, s, r, k, z;
	for (int i=1; i<=number_of_unknowns; i++)
		x(i)=laplace->u(matrix_index_to_cell_index(i));
	laplace->Find_Tolerance(b); // needs to happen after b is completely set up
	DOMAIN_ITERATOR_THREADED_ALPHA<PCG_SPARSE_THREADED<TV>,TV> threaded_iterator(
			laplace->grid.Domain_Indices(1), laplace->thread_queue, 1, 1, 2, 1);
	
	// Initialize the projection driver with MPI. 
	// [TODO] Remove MPI!
	MPI::Intracomm* comm;
	comm=&(*laplace->laplace_mpi->communicators)(color);
	pcg_mpi	= new NIMBUS_PCG_SPARSE_MPI<GRID<TV> >(laplace->laplace_mpi->local_pcg,*comm,laplace->laplace_mpi->partitions(color));

	// pcg_mpi.Parallel_Solve(A,x,b,laplace->tolerance);
	projection_data->tolerance = laplace->tolerance;	
}

template<class TV> void PROJECTION_DRIVER<TV>::WindUpForOneRegion() {
	int color = 1;
	LAPLACE_UNIFORM<GRID<TV> >* laplace = example.projection.elliptic_solver;
	ARRAY<TV_INT> &matrix_index_to_cell_index =
			(*(projection_data->matrix_index_to_cell_index_array))(color);
	VECTOR_ND<T> &x = *(projection_data->x);
	int number_of_unknowns=matrix_index_to_cell_index.m;
	for (int i=1; i<=number_of_unknowns; i++) {
		TV_INT cell_index=matrix_index_to_cell_index(i);
		laplace->u(cell_index)=x(i);
	}
}

template<class TV> void PROJECTION_DRIVER<TV>::Write_Substep(
		const std::string& title, const int substep, const int level) {
	if (level<=example.write_substeps_level) {
		example.frame_title=title;
		std::stringstream ss;
		ss<<"Writing substep ["<<example.frame_title<<"]: output_number="
				<<output_number+1<<", time="<<time<<", frame="<<current_frame
				<<", substep="<<substep<<std::endl;
		LOG::filecout(ss.str());
		Write_Output_Files(++output_number);
		example.frame_title="";
	}
}

template<class TV> void PROJECTION_DRIVER<TV>::Write_Output_Files(
		const int frame) {
	FILE_UTILITIES::Create_Directory(example.output_directory);
	FILE_UTILITIES::Create_Directory(example.output_directory
			+STRING_UTILITIES::string_sprintf("/%d", frame));
	FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
	FILE_UTILITIES::Write_To_Text_File(example.output_directory
			+STRING_UTILITIES::string_sprintf("/%d/frame_title", frame),
			example.frame_title);
	if (frame==example.first_frame)
		FILE_UTILITIES::Write_To_Text_File(example.output_directory
				+"/common/first_frame", frame, "\n");
	example.Write_Output_Files(frame);
	FILE_UTILITIES::Write_To_Text_File(example.output_directory
			+"/common/last_frame", frame, "\n");
}

template class PROJECTION_DRIVER<VECTOR<float,2> >;
template class PROJECTION_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_DRIVER<VECTOR<double,2> >;
template class PROJECTION_DRIVER<VECTOR<double,3> >;
#endif
