//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PRECONDITIONER
//##################################################################### 
#ifndef __CONTACT_PRECONDITIONER__
#define __CONTACT_PRECONDITIONER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PROJECTED_GAUSS_SEIDEL.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>

namespace PhysBAM{

namespace CONTACT_PRECONDITIONER
{
// v1 = v0 + M^-1 C^t lambda
// C v1 >= b
// C M^-1 C^t lambda >= b - C v1

template<class TV>
class PRECONDITIONER
{
public:
    typedef typename TV::SCALAR T;
    typedef CONTACT<TV> T_CONTACT;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
    static const int d=TV::dimension;
    static const int d_spin=TV::SPIN::dimension;
    static const int d_rigid=d+d_spin;

    ARRAY<int> body_to_dynamic;

    SPARSE_MATRIX_FLAT_MXN<T> A; //C M^-1 C^t //ISN'T THIS SQUARE?
    VECTOR_ND<T> a; //b - C v0
    VECTOR_ND<T> lambda;

    void Build_Preconditioner(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<T_CONTACT>& contacts)
    {
        LOG::SCOPE scope("PRECONDITIONER::Build_Preconditioner");

        int n_bodies=rigid_body_collection.rigid_body_particle.array_collection->Size();
        int n_contacts=contacts.m;

        {std::stringstream ss;ss << "bodies " << n_bodies << " contacts " << n_contacts << std::endl;LOG::filecout(ss.str());}

        int n_bodies_dynamic=0;
        body_to_dynamic.Resize(n_bodies);
        for(int i=1;i<=n_bodies;i++)
            if(rigid_body_collection.Is_Active(i)){
                if(rigid_body_collection.Rigid_Body(i).Has_Infinite_Inertia()) body_to_dynamic(i)=0;
                else body_to_dynamic(i)=++n_bodies_dynamic;}
            else body_to_dynamic(i)=-1;

        SPARSE_MATRIX_FLAT_MXN<T> Mi; //mass inverse
        SPARSE_MATRIX_FLAT_MXN<T> C,Ct;

        ARRAY<int> Mi_row_lengths(d_rigid*n_bodies_dynamic);
        for(int i=1;i<=n_bodies;i++){int index=body_to_dynamic(i);
            if(index>0){
                for(int j=1;j<=d;j++) Mi_row_lengths((index-1)*d_rigid+j)=1;
                for(int j=1;j<=d_spin;j++) Mi_row_lengths((index-1)*d_rigid+j+d)=3;}}
        
        Mi.Set_Row_Lengths(Mi_row_lengths);
        Mi.n=d_rigid*n_bodies_dynamic;
        for(int i=1;i<=n_bodies;i++){int index=body_to_dynamic(i);
            if(index>0){
                RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
                T mass_inverse=1.0/body.Mass();
                T_WORLD_SPACE_INERTIA_TENSOR inertia_tensor_inverse=body.World_Space_Inertia_Tensor_Inverse();
                int base=(index-1)*d_rigid;
                for(int j=1;j<=d;j++) Mi.Set_Element(base+j,base+j,mass_inverse);
                for(int j=1;j<=d_spin;j++) for(int k=1;k<=d_spin;k++) Mi.Set_Element(base+d+j,base+d+k,inertia_tensor_inverse(j,k));}}
        
        ARRAY<int> C_row_lengths(n_contacts);
        for(int i=1;i<=n_contacts;i++) C_row_lengths(i)=(((body_to_dynamic(contacts(i).id(1))>0)?1:0)+(body_to_dynamic(contacts(i).id(2))>0?1:0))*d_rigid;
        C.Set_Row_Lengths(C_row_lengths);
        C.n=d_rigid*n_bodies_dynamic;
        for(int i=1;i<=n_contacts;i++){
            T_CONTACT& contact=contacts(i);
            for(int j=1;j<=2;j++){int index=body_to_dynamic(contact.id(j));
                if(index){int base=(index-1)*d_rigid;
                    for(int k=1;k<=d;k++) C.Set_Element(i,base+k,contact.normal_constraint(j).linear(k));
                    for(int k=1;k<=d_spin;k++) C.Set_Element(i,base+k+d,contact.normal_constraint(j).angular(k));}}}

        C.Transpose(Ct);
        A=C*Mi*Ct;
        
        a.Resize(n_contacts);
        for(int i=1;i<=n_contacts;i++){
            T_CONTACT& contact=contacts(i);
            T b=contact.normal_relative_velocity;
            for(int j=1;j<=2;j++) b-=TWIST<TV>::Dot_Product(contact.normal_constraint(j),rigid_body_collection.Rigid_Body(contact.id(j)).Twist());
            a(i)=b;}

        lambda.Resize(n_contacts);
    }

    void Update_Velocities(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<T_CONTACT>& contacts)
    {
        LOG::SCOPE scope("PRECONDITIONER::Update_Velocities");

        int n_bodies=rigid_body_collection.rigid_body_particle.array_collection->Size();
        int n_contacts=contacts.m;
        
        for(int i=1;i<=n_contacts;i++){
            T_CONTACT& contact=contacts(i);
            for(int j=1;j<=2;j++){
                RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(contact.id(j));
                body.Twist()+=contact.inverse_mass_times_normal_constraint(j)*lambda(i);}}
        
        for(int i=1;i<=n_bodies;i++)
            if(rigid_body_collection.Is_Active(i) && !rigid_body_collection.Rigid_Body(i).Has_Infinite_Inertia())
                rigid_body_collection.Rigid_Body(i).Update_Angular_Momentum();
    }

    T Residual(SPARSE_MATRIX_FLAT_MXN<T>& A,VECTOR_ND<T>& a,VECTOR_ND<T>& lambda)
    {
        T residual_max=0;
        int n_contacts=A.m;
        VECTOR_ND<T> r(n_contacts);
        A.Times(lambda,r);
        r-=a;
        for(int i=1;i<=n_contacts;i++){
            T residual=(r(i)<0)?(-r(i)):fabs(r(i)*lambda(i));
            if(residual>residual_max)
                residual_max=residual;}
        return residual_max;
    }

    void Solve(T tolerance,int iteration_maximum,bool recursive=true)
    {
        LOG::SCOPE scope("PRECONDITIONER::Solve");

        int iteration;
        for(iteration=0;Residual(A,a,lambda)>tolerance && iteration<iteration_maximum;iteration++)
            Solve(A,a,lambda,tolerance,iteration_maximum);

        {std::stringstream ss;ss << "outer iterations " << iteration << std::endl;LOG::filecout(ss.str());}
    }

    void Sparsify(SPARSE_MATRIX_FLAT_MXN<T>& A,SPARSE_MATRIX_FLAT_MXN<T>& A_coarse,T tolerance)
    {
        int n_contacts=A.m;

        VECTOR_ND<T> diagonal(n_contacts);
        for(int i=1;i<=n_contacts;i++) for(int j=A.offsets(i);j<A.offsets(i+1);j++) if(A.A(j).j==i){diagonal(i)=A.A(j).a;break;}

        VECTOR_ND<T> diagonal_offset(n_contacts);
        diagonal_offset.Fill(0);
        ARRAY<int> A_coarse_row_lengths(n_contacts);
        for(int i=1;i<=n_contacts;i++){
            int c=0;
            for(int j=A.offsets(i);j<A.offsets(i+1);j++){
                int index=A.A(j).j;
                if(index==i || fabs(A.A(j).a)/sqrt(diagonal(i)*diagonal(index))>tolerance) c++;
                else diagonal_offset(i)+=fabs(A.A(j).a);}
            A_coarse_row_lengths(i)=c;}
        A_coarse.Set_Row_Lengths(A_coarse_row_lengths);
        A_coarse.n=n_contacts;
        for(int i=1;i<=n_contacts;i++){
            for(int j=A.offsets(i);j<A.offsets(i+1);j++){
                int index=A.A(j).j;
                //{std::stringstream ss;ss << "setting value " << i << " " << index << " " << A.A(j).a << std::endl;LOG::filecout(ss.str());}
                if(index==i || fabs(A.A(j).a)/sqrt(diagonal(i)*diagonal(index))>tolerance) A_coarse.Set_Element(i,index,A.A(j).a);}
            A_coarse.Add_Element(i,i,diagonal_offset(i));}

        OCTAVE_OUTPUT<typename TV::SCALAR> output("matrices.m");
        output.Write("A",A);
        output.Write("A_coarse",A_coarse);
        
        //{std::stringstream ss;ss << "A = [ " << A << " ]" << std::endl;LOG::filecout(ss.str());}
        //{std::stringstream ss;ss << "A_coarse = [ " << A_coarse << " ]" << std::endl;LOG::filecout(ss.str());}
    }

    void Build_Coarse_System(SPARSE_MATRIX_FLAT_MXN<T>& A,VECTOR_ND<T>& a,SPARSE_MATRIX_FLAT_MXN<T>& A_coarse,VECTOR_ND<T>& a_coarse)
    {
        Sparsify(A,A_coarse,0.05);
        a_coarse=a;
    }

    void Solve(SPARSE_MATRIX_FLAT_MXN<T>& A,VECTOR_ND<T>& a,VECTOR_ND<T>& x,T tolerance,int iteration_maximum)
    {
        int n_contacts=A.m;

        SPARSE_MATRIX_FLAT_MXN<T> A_coarse;
        VECTOR_ND<T> a_coarse;
        Build_Coarse_System(A,a,A_coarse,a_coarse);

        {std::stringstream ss;ss << "residual 0 " << Residual(A,a,lambda) << std::endl;LOG::filecout(ss.str());}

        PROJECTED_GAUSS_SEIDEL::Solve(A,a,lambda,tolerance,10);

        VECTOR_ND<T> r(n_contacts);
        A.Times(lambda,r);
        r=a_coarse-r;
        VECTOR_ND<T> e(n_contacts);

        {std::stringstream ss;ss << "residual 1 " << Residual(A,a,lambda) << std::endl;LOG::filecout(ss.str());}

        PROJECTED_GAUSS_SEIDEL::Solve(A_coarse,r,e,tolerance,100);
        lambda+=e;

        {std::stringstream ss;ss << "residual 2 " << Residual(A,a,lambda) << std::endl;LOG::filecout(ss.str());}

        PROJECTED_GAUSS_SEIDEL::Solve(A,a,lambda,tolerance,10);

        {std::stringstream ss;ss << "residual 3 " << Residual(A,a,lambda) << std::endl;LOG::filecout(ss.str());}

        //exit(0);
    }
};

}
}

#endif
