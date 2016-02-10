//#####################################################################
// Copyright 2011, Valeria Nikolaenko
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_ANIMATED_VISUALIZATION_INTERFACE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_ANIMATED_VISUALIZATION_INTERFACE__
#define __OPTIX_ANIMATED_VISUALIZATION_INTERFACE__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/SELECTION_RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_RIGID_GEOMETRY_COLLECTION.h>
#include <climits>

namespace PhysBAM {

template<class TV>
class OPTIX_ANIMATED_VISUALIZATION_INTERFACE {
    typedef typename TV::SCALAR TS;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
public:
    virtual void Set_Frame(int frame_input) = 0;
    virtual bool Has_Current_Selection() { return false; }
    virtual SELECTION_RIGID_BODY<TV>* Get_Current_Selection() {return NULL;}

    /*virtual void Set_MAC_Velocities_Simulated(ARRAY<T,FACE_INDEX<TV::dimension> >* face_velocities_simulated);
    virtual void Set_Scalar_Values_Simulated(ARRAY<T,TV_INT>* densities_simulated);
    virtual void Set_Rigid_Bodies_Simulated(RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection_simulated);

    virtual ARRAY<T,FACE_INDEX<TV::dimension> >* Get_MAC_Velocities_Simulated();

    virtual RIGID_GEOMETRY_COLLECTION<TV>* Get_Rigid_Bodies_Simulated();
    */
    // virtual SCALAR_COMPONENT<TV>* Get_Scalar_Component() = 0;
    virtual void Set_Grid(GRID<TV> _grid) = 0;
    virtual GRID<TV>& Get_Grid() = 0;

    virtual void Set_Scalar_Values_Simulated(ARRAY<TS,TV_INT>* densities_simulated) = 0;
    virtual ARRAY<TS,TV_INT>* Get_Scalar_Values_Simulated() = 0;

    virtual void Set_Scalar_Values_Upsample_Scale(int upsample_scale) = 0;
    virtual int Get_Upsample_Scale() = 0;

    virtual void Set_Rigid_Bodies_Simulated(RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection_simulated) = 0;
    virtual RIGID_GEOMETRY_COLLECTION<TV>* Get_Rigid_Bodies_Simulated() = 0;

    virtual void Rigid_Geometry_Updated() = 0;
};

}

#endif
#endif
