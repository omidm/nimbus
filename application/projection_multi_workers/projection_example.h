/*
 * Minor change from original PhysBAM code "PROJECTION_EXAMPLE.h".
 * I think these data structures should have relatively small size in future,
 * only containing metadata such as time step and grid size.
 * It is expected to be in the profile data abstraction.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#ifndef __PROJECTION_EXAMPLE__
#define __PROJECTION_EXAMPLE__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

template<class TV> class PROJECTION_EXAMPLE {
	typedef typename TV::SCALAR T;
	typedef typename TV::template REBIND<int>::TYPE TV_INT;
	enum workaround1 {d=TV::m};
public:
	STREAM_TYPE stream_type;
	T initial_time;
	int first_frame;
	T frame_rate;
	int restart;
	std::string frame_title;
	int write_substeps_level;
	bool write_debug_data;
	std::string output_directory;

	GRID<TV> mac_grid;
	MPI_UNIFORM_GRID<GRID<TV> > *mpi_grid;
	THREAD_QUEUE* thread_queue;
	PROJECTION_UNIFORM<GRID<TV> > projection;
	ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
	BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
	BOUNDARY_UNIFORM<GRID<TV>,T> *boundary;
	VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;

	PROJECTION_EXAMPLE(const STREAM_TYPE stream_type_input,
			int number_of_threads=1);
	virtual ~PROJECTION_EXAMPLE();

	T Time_At_Frame(const int frame) const {
		return initial_time+(frame-first_frame)/frame_rate;
	}

	void Initialize_Grid(TV_INT counts, RANGE<TV> domain) {
		mac_grid.Initialize(counts, domain, true);
	}

	void Initialize_Fields() {
		for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid); iterator.Valid(); iterator.Next())
			face_velocities(iterator.Full_Index())
					=(T)iterator.Face_Index()(iterator.Axis())
							/(T)mac_grid.counts(iterator.Axis());

		/*
		for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid); iterator.Valid(); iterator.Next()) {
			face_velocities(iterator.Full_Index())
					=(T)iterator.Face_Index()(iterator.Axis())
							/(T)mac_grid.counts(iterator.Axis())
							+ ((iterator.Axis()==2) ? (T)mpi_grid->rank/2 : 0);
			face_velocities(iterator.Full_Index())
					*= face_velocities(iterator.Full_Index());
		}
		*/
	}

	virtual void Write_Output_Files(const int frame);
	virtual void Read_Output_Files(const int frame);
	virtual void Set_Boundary_Conditions(const T time);
};
} // namespace PhysBAM

#endif
