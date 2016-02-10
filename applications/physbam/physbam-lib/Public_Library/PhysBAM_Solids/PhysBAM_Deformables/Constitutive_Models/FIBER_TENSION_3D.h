//#####################################################################
// Copyright 2005, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FIBER_TENSION_3D
//##################################################################### 
#ifndef __FIBER_TENSION_3D__
#define __FIBER_TENSION_3D__

#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
namespace PhysBAM{

template<class T>
class FIBER_TENSION_3D:public ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ANISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;

    T fiber_p1,fiber_p2,fiber_p3,fiber_p4,fiber_cutoff;
    ARRAY<T> fiber_max_stress;
    T failure_threshold;
    ARRAY<T>& muscle_activations;
    ARRAY<ARRAY<int> > tet_muscles;
    ARRAY<ARRAY<VECTOR<T,3> > > tet_fibers;
    ARRAY<ARRAY<T> > tet_densities;
    ARRAY<ARRAY<T> > tension;
    bool use_force_length;

    FIBER_TENSION_3D(ARRAY<T>& muscle_activations_input,const ARRAY<T>* peak_isometric_stress,const T failure_threshold_input,const bool use_force_length_input)
        :failure_threshold(failure_threshold_input),muscle_activations(muscle_activations_input),use_force_length(use_force_length_input)
    {
        fiber_p1=(T).05;
        fiber_p2=(T)6.6;
        fiber_cutoff=(T)1.4;
        T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
        if(peak_isometric_stress) {assert(peak_isometric_stress->m==muscle_activations.m);fiber_max_stress=*peak_isometric_stress;} 
        else {fiber_max_stress.Resize(muscle_activations.m);ARRAY<T>::Copy((T)3e5,fiber_max_stress);}
        fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);
        fiber_p4=fiber_p1*(exp(cutoff_scaled)*(1-fiber_p2*fiber_cutoff)+fiber_p2-1);
    }

    void Initialize_Fiber_Data(const STRAIN_MEASURE<TV,3>& strain_measure,const ARRAY<ARRAY<int> >& muscle_tets,const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,
        const ARRAY<ARRAY<T> >& muscle_densities)
    {int n=strain_measure.tetrahedron_mesh.tetrahedrons.m;
    tet_muscles.Resize(n);tet_fibers.Resize(n);tension.Resize(n);tet_densities.Resize(n);
    for(int m=1;m<=muscle_tets.m;m++) for(int t=1;t<=muscle_tets(m).m;t++){
        tet_muscles(muscle_tets(m)(t)).Append(m);tet_densities(muscle_tets(m)(t)).Append(muscle_densities(m)(t));
        tet_fibers(muscle_tets(m)(t)).Append(strain_measure.F(muscle_tets(m)(t)).Transpose_Times(muscle_fibers(m)(t)));
        tension(muscle_tets(m)(t)).Append(0);}}

    T Tension(const int tet_index,const int tet_muscle_index,const T stretch) const
    {if(use_force_length){
        T strain=stretch-1,strain_abs=abs(strain),activation=muscle_activations(tet_muscles(tet_index)(tet_muscle_index)),density=tet_densities(tet_index)(tet_muscle_index);
        T active_tension=0,passive_tension=0,scale=(T)25/(T)6;
        if(stretch>fiber_cutoff)passive_tension=fiber_p3*stretch+fiber_p4;else if(stretch>1)passive_tension=fiber_p1*(exp(fiber_p2*strain)-fiber_p2*strain-1);
        if(strain_abs<.4)active_tension=activation*density*(1-scale*sqr(strain));else if(strain_abs<.6)active_tension=2*scale*activation*density*sqr(strain_abs-(T).6);
        return fiber_max_stress(tet_muscles(tet_index)(tet_muscle_index))*(active_tension+passive_tension);}
    else
        return muscle_activations(tet_muscles(tet_index)(tet_muscle_index))*tet_densities(tet_index)(tet_muscle_index)*fiber_max_stress(tet_muscles(tet_index)(tet_muscle_index));}
        
    MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const T scale,const int tetrahedron_index) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold);MATRIX<T,3> fiber_stress;
    for(int m=1;m<=tet_muscles(tetrahedron_index).m;m++){
        VECTOR<T,3> fiber=V.Transpose_Times(tet_fibers(tetrahedron_index)(m)),F_fiber=F_threshold*fiber;
        fiber_stress+=(tension(tetrahedron_index)(m)/F_fiber.Magnitude())*MATRIX<T,3>::Outer_Product(F_fiber,fiber);}
    return scale*fiber_stress;}

    void Update_State_Dependent_Auxiliary_Variables(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const int tetrahedron_index) PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold);
    for(int m=1;m<=tet_muscles(tetrahedron_index).m;m++){
        VECTOR<T,3> fiber=V.Transpose_Times(tet_fibers(tetrahedron_index)(m)),F_fiber=F_threshold*fiber;T stretch=F_fiber.Magnitude();
        tension(tetrahedron_index)(m)=Tension(tetrahedron_index,m,stretch);}}

//#####################################################################
};

template<class T> FINITE_VOLUME<VECTOR<T,3>,3>*
Create_Fiber_Tension(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<ARRAY<int> >& muscle_tets,
    const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,const ARRAY<ARRAY<T> >& muscle_densities,ARRAY<T>& muscle_activations,const bool use_force_length=true,
    const ARRAY<T>* peak_isometric_stress=0,const T failure_threshold=(T).25,const bool verbose=false)
{
    FIBER_TENSION_3D<T>* constitutive_model=new FIBER_TENSION_3D<T>(muscle_activations,peak_isometric_stress,failure_threshold,use_force_length);
    FINITE_VOLUME<VECTOR<T,3>,3>* fvm=Create_Finite_Volume(tetrahedralized_volume,constitutive_model,false,FLT_MAX,false,verbose,0);
    fvm->use_velocity_dependent_forces=false;
    constitutive_model->Initialize_Fiber_Data(fvm->strain_measure,muscle_tets,muscle_fibers,muscle_densities);
    return fvm;
}

}
#endif
