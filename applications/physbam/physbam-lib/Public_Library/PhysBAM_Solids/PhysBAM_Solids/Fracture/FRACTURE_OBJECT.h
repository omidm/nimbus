//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ronald Fedkiw, Frank Losasso, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_OBJECT
//#####################################################################
#ifndef __FRACTURE_OBJECT__
#define __FRACTURE_OBJECT__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV,int d> class LEVELSET_GRAIN_BOUNDARIES;

template<class TV,int d>
class FRACTURE_OBJECT:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    typedef typename MATRIX_POLICY<TV>::DIAGONAL_MATRIX T_DIAGONAL_MATRIX;
    typedef typename MESH_POLICY<d>::MESH T_MESH;
public:
    T eigenvector_coefficient; // post processing fracture normal
    T fracture_bias_direction_coefficient; // post processing fracture normal
    T fracture_bias_propagation_coefficient; // post processing fracture normal
    T fracture_bias_propagation; // control
    ARRAY<T> fracture_bias_magnitude; // control
    ARRAY<VECTOR<T,d> > fracture_bias_direction; // control
    ARRAY<T> fracture_bias_stress_scaling; // control
    ARRAY<int> fracture_phi_index;
    ARRAY<VECTOR<T,d+1> > phi; // lower dimensional    
    int number_of_fracture_initiations;
    T_EMBEDDED_OBJECT& embedded_object;
private:
    PARTICLES<TV> reference_particles;
    T_MESH reference_mesh;
public:
    T_SIMPLICIAL_OBJECT reference_simplicial_object;
    VECTOR<T,d> fracture_threshold;
    VECTOR<T,d> compressive_threshold;
    int max_number_of_cuts;
    int number_of_smoothing_passes;
    int number_of_fracture_initiation_points;
    bool bias_stress;
    T fracture_quality_threshold;
    T extra_levelset_propagation_direction_scaling;
    bool force_edge_connected_fracture;
    ARRAY<int> corresponding_node_in_reference;
    ARRAY<int> corresponding_simplex_in_reference;
protected:
    RANDOM_NUMBERS<T> random_numbers; // control
public:

    // Optional Parameters
    ARRAY<TV>* initiation_point_positions; // reference positions of initiations points.
    ARRAY<TV>* initiation_point_reference_seed_positions;   // hard code initiation point positions
    ARRAY<T>* initiation_point_radii; // radii of initiation points.

    FRACTURE_OBJECT(T_EMBEDDED_OBJECT& embedded_object_input);
    virtual ~FRACTURE_OBJECT();

    void Set_Fracture_Bias_Magnitude(const T fracture_bias_magnitude_input=0)
    {fracture_bias_magnitude.Resize(reference_mesh.elements.m,false);ARRAYS_COMPUTATIONS::Fill(fracture_bias_magnitude,fracture_bias_magnitude_input);}

    void Set_Fracture_Bias_Stress_Scaling(const T stress_scale_input=1)
    {fracture_bias_stress_scaling.Resize(reference_mesh.elements.m,false);ARRAYS_COMPUTATIONS::Fill(fracture_bias_stress_scaling,stress_scale_input);}

    void Set_Fracture_Bias_Propagation(const T fracture_bias_propagation_input=0)
    {fracture_bias_propagation=fracture_bias_propagation_input;}

    int Fracture_Phi_Index(const int t) const
    {return fracture_phi_index(corresponding_simplex_in_reference(t));}

    void Set_Phi(const int t,const VECTOR<T,d+1>& phi_element)
    {assert(fracture_phi_index.m==reference_mesh.elements.m);
    int ref_t=corresponding_simplex_in_reference(t);
    int index=fracture_phi_index(ref_t);if(index) phi(index)=phi_element;else fracture_phi_index(ref_t)=phi.Append(phi_element);}

    VECTOR<T,d+1> Get_Phi(const int t) const
    {return phi(Fracture_Phi_Index(t));}

    int Positive_Count(const int t) const
    {VECTOR<T,d+1> phi_element=Get_Phi(t);
    int count=0;for(int i=1;i<=d+1;i++) if(phi_element[i]>0) count++;return count;}

    T Phi_In_Simplex(const int node,const int simplex) const
    {assert(embedded_object.simplicial_object.mesh.Node_In_Simplex(node,simplex));
    return Get_Phi(simplex)[embedded_object.simplicial_object.mesh.elements(simplex).Find(node)];}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,number_of_fracture_initiations,corresponding_node_in_reference,corresponding_simplex_in_reference);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,number_of_fracture_initiations,corresponding_node_in_reference,corresponding_simplex_in_reference);}

//#####################################################################
    virtual void Rebuild_Embedded_Object(ARRAY<int>& map_to_old_particles,ARRAY<int>& map_to_old_embedded_particles,ARRAY<int>& map_to_old_simplices,const bool verbose=true);
    int Fracture_Where_High_Stress(ARRAY<T_SYMMETRIC_MATRIX>& sigma,ARRAY<TV>& spatial_fracture_bias_direction);
    TV Rest_Position_Of_Material_Surface_Particle(const int boundary_particle);
    void Set_Random_Fracture_Bias_Stress_Scaling_Constant(const T fracture_threshold,const int averaging_iterations); // TODO: this function is apparently deprecated
    virtual void Add_Cut_Based_On_Phi(const int element,const VECTOR<T,d+1>& tetrahedron_phi){};
private:
    virtual void Add_Cut(const int element,const TV& fracture_normal)=0;
protected:
    bool Initiation_Point(const int element);
//#####################################################################
};
}
#include <PhysBAM_Solids/PhysBAM_Solids/Read_Write/Fracture/READ_WRITE_FRACTURE_OBJECT.h>
#endif
