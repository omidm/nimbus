//#####################################################################
// Copyright 2004-2005, Igor Neverov, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_CONTROL_SET
//#####################################################################
#ifndef __FACE_CONTROL_SET__    
#define __FACE_CONTROL_SET__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;
template<class TV> class TWIST;

template<class T>
class FACE_CONTROL_SET
{
    typedef VECTOR<T,3> TV;
public:
    enum TYPE {ACTIVATION,ATTACHMENT_FRAME};
    FACE_CONTROL_SET()
    {}

    virtual ~FACE_CONTROL_SET()
    {}
    
//#####################################################################
    virtual int Size() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TYPE Type() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T operator()(const int control_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T& operator()(const int control_id) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Control_Active(const int control_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool& Control_Active(const int control_id) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Positions_Determined_Kinematically(const int control_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Force_Derivative(ARRAY<TV>& dFdl,ARRAY<TWIST<TV> >& dFrdl,const int control_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Position_Derivative(ARRAY<TV>& dXdl,const int control_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Set_Attachment_Positions(ARRAY<TV>&X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Save_Controls() {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Kinematically_Update_Positions(ARRAY<TV>&X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Kinematically_Update_Jacobian(ARRAY<TV>&dX) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Penalty() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Penalty_Gradient(const int control_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Penalty_Hessian(const int control_id1,const int control_id2) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Print_Diagnostics(std::ostream& output) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Project_Parameters_To_Allowable_Range(const bool active_controls_only=false) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Interpolate(const T interpolation_fraction) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Scale(const T scale) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Distance(const ARRAY<T>& weights,int base) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Identity(const int control_id) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Maximal_Controls() {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
