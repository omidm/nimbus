//#####################################################################
// Copyright 2002-2012, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Igor Neverov, Craig Schroeder, Rahul Sheth, Tamar Shinar, Eftychios Sifakis, Joseph Teran, Rachel Weinstein, Zahid Hossain
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_ND
//#####################################################################
#ifndef __VECTOR_ND__
#define __VECTOR_ND__

#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_EXPRESSION.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
namespace PhysBAM{

template<class T> struct IS_ARRAY_VIEW<VECTOR_ND<T> > {static const bool value=false;};

template<class T>
class VECTOR_ND:public VECTOR_BASE<T,VECTOR_ND<T> >
{
public:
    typedef T SCALAR;typedef VECTOR_BASE<T,VECTOR_ND<T> > BASE;
    template<class T2> struct REBIND{typedef VECTOR_ND<T2> TYPE;};
    typedef T ELEMENT;
    typedef int INDEX;

    int n; // size of the n vector
    T* x; // pointer to the n vector
private:
    bool owns_data; // whether or not we own x
public:
    size_t dynamic_size() {
      return sizeof(T) * n;
    }
    size_t memory_size() {
      return sizeof(*this) + dynamic_size();
    }

    explicit VECTOR_ND(const int n_input=0,const bool initialize=true)
        :n(n_input),owns_data(true)
    {
        x=new T[n];if(initialize) for(int k=0;k<n;k++) x[k]=0;
    }

    explicit VECTOR_ND(const INITIAL_SIZE n_input)
        :n(Value(n_input)),owns_data(true)
    {
        x=new T[n];for(int k=0;k<n;k++) x[k]=0;
    }

    VECTOR_ND(const VECTOR_ND& v)
        :n(v.n),owns_data(true)
    {
        x=new T[n];for(int k=0;k<n;k++) x[k]=v.x[k];
    }

    template<class T_VECTOR>
    explicit VECTOR_ND(const VECTOR_BASE<T,T_VECTOR>& vector)
        :owns_data(true)
    {
        const T_VECTOR& derived=vector.Derived();
        n=derived.Size();x=new T[n];
        for(int k=1;k<=n;k++) (*this)(k)=derived(k);
    }

    template<class T2,class T_VECTOR>
    explicit VECTOR_ND(const VECTOR_BASE<T2,T_VECTOR>& vector)
        :owns_data(true)
    {
        const T_VECTOR& derived=vector.Derived();
        n=derived.Size();x=new T[n];
        for(int k=1;k<=n;k++) (*this)(k)=(T)derived(k);
    }

    ~VECTOR_ND()
    {if(owns_data) delete[] x;}

    void Set_Subvector_View(const VECTOR_ND& v,const INTERVAL<int>& range)
    {assert(1<=range.min_corner && range.min_corner<=range.max_corner+1 && range.max_corner<=v.n);
    if(owns_data) delete[] x;owns_data=false;n=range.max_corner-range.min_corner+1;x=v.x+range.min_corner-1;}

    VECTOR_ND& operator=(const VECTOR_ND& source)
    {if(n!=source.n){assert(owns_data);n=source.n;delete[] x;x=new T[n];}
    for(int k=0;k<n;k++) x[k]=source.x[k];return *this;}

    template<class T_VECTOR>
    VECTOR_ND& operator=(const VECTOR_BASE<T,T_VECTOR>& vector)
    {const T_VECTOR& derived=vector.Derived();int size=derived.Size();
    if(n!=size){assert(owns_data);n=size;delete[] x;x=new T[n];}
    for(int k=1;k<=n;k++) (*this)(k)=derived(k);return *this;}

    int Size() const
    {return n;}

    T& operator()(const int i)
    {assert(i>=1 && i<=n);return x[i-1];}

    const T& operator()(const int i) const
    {assert(i>=1 && i<=n);return x[i-1];}

    T& operator[](const int i)
    {assert(i>=1 && i<=n);return x[i-1];}

    const T& operator[](const int i) const
    {assert(i>=1 && i<=n);return x[i-1];}


    T* Get_Array_Pointer()
    {return x;}

    const T* Get_Array_Pointer() const
    {return x;}

    void Resize(const int n_new,const bool initialize=true)
    {assert(owns_data);
    if(n_new==n)return;
    T* x_new=new T[n_new];
    int n1=PhysBAM::min(n,n_new);for(int i=0;i<n1;i++) x_new[i]=x[i];if(initialize) for(int i=n1;i<n_new;i++) x_new[i]=0;
    delete[] x;x=x_new;n=n_new;}

    bool Owns_Data() const
    {return owns_data;}

    void Sort_Descending(){std::sort(x,x+n,std::greater<T>());}
    void Sort_Ascending(){std::sort(x,x+n,std::less<T>());}

    int Arg_Max() const
    {assert(1 <= n);
    T max = x[0];int arg_max = 1;
    for(int i=0;i<n;++i) if( max < x[i] ) { max = x[i]; arg_max = i+1; };
    return arg_max;}

    VECTOR_ND<T> operator-( const VECTOR_ND<T> &v) const {
        assert( v.Size() == Size() );
        VECTOR_ND<T> result( Size() );
        for(int i=1;i<=Size();i++) result(i) = (*this)(i) - v(i);
        return result;
    }

    bool All_Less_Equal( const VECTOR_ND<T> &v) const {
        for(int i=1;i<Size();i++) { if( (*this)(i) > v(i) ) return false; }
        return true;
    }
//#####################################################################
};
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#endif
#endif
