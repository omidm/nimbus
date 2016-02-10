//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_COLLISIONS
//#####################################################################
#ifndef __BW_COLLISIONS__
#define __BW_COLLISIONS__    

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>

namespace PhysBAM{

template<class T1,class T2> class PAIR;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class BW_BACKWARD_EULER_SYSTEM;
template<class T> class KRYLOV_VECTOR_BASE;

template<class TV>
class BW_COLLISIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef GENERALIZED_VELOCITY<TV> VECTOR_T;
public:
    ARRAY<PAIR<int,int> > cloth_body_constraints;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;

    BW_COLLISIONS(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input);
    virtual ~BW_COLLISIONS();

//#####################################################################
    void Detect_Cloth_Body_Contact();
    void Remove_Separating_Cloth_Body_Contacts(BW_BACKWARD_EULER_SYSTEM<TV>& system,KRYLOV_VECTOR_BASE<T>& R,KRYLOV_VECTOR_BASE<T>& B,KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& Q);
//#####################################################################
};
}
#endif
