//#####################################################################
// Copyright 2012, Rahul Sheth, Zahid Hossain
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STATE_VECTOR
//#####################################################################
// This class helps map model reduction state vectors to PhysBAM
// data structures.
//#####################################################################
#ifndef __STATE_VECTOR__
#define __STATE_VECTOR__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>

namespace PhysBAM{

template<class TV>
class STATE_VECTOR
{
public:
    typedef typename TV::SCALAR T;
    typedef T SCALAR;

    int length; //total size of this vector
    int nodes; // number of nodes we have the state for
    
    int dimension; //1, 2. or 3

    T** data;

    STATE_VECTOR(GEOMETRY_PARTICLES<TV>& geometry_particles_in)
    {
        dimension = TV::dimension;
        int doubleDimension = dimension*2;

        length=geometry_particles_in.X.m*doubleDimension;
        nodes=geometry_particles_in.X.m;

        data=new T*[length];
        for(int i=1;i<=length;i+=doubleDimension)
        {
            //TODO: There's probably a better way to do this...
            data[i-1]=&geometry_particles_in.X((i/(doubleDimension))+1)(1);
            data[i]=&geometry_particles_in.X((i/(doubleDimension))+1)(2);
            if(dimension==2)
            {
                data[i+1]=&geometry_particles_in.V((i/doubleDimension)+1)(1);
                data[i+2]=&geometry_particles_in.V((i/doubleDimension)+1)(2);
            }
            else if (dimension==3)
            {
                data[i+1]=&geometry_particles_in.X((i/doubleDimension)+1)(3);
                data[i+2]=&geometry_particles_in.V((i/doubleDimension)+1)(1);
                data[i+3]=&geometry_particles_in.V((i/doubleDimension)+1)(2);
                data[i+4]=&geometry_particles_in.V((i/doubleDimension)+1)(3);
            }
            else 
            {
                PHYSBAM_FATAL_ERROR("STATE_VECTOR: Does not support dimension!");
            }
        }
    }

    STATE_VECTOR(VECTOR_ND<T>& input_vector)
    {
        dimension=TV::dimension;

        length=input_vector.Size();
        nodes=input_vector.Size()/dimension;

        data=new T*[length];
        for(int i=1;i<=length;i++)
        {
            data[i-1]=&input_vector(i);
        }
    }

    ~STATE_VECTOR()
    {
        delete data;
    }

    inline int Total_Size() const {return length;}
    
    inline int Number_Of_Nodes() const {return nodes;}
    
    inline int Size() const {return nodes;}

    T& operator()(const int i)
    {
        assert(1<=i && i<=length);
        return *data[i-1];
    }

    T& operator()(const int i) const
    {
        assert(1<=i && i<=length);
        return *data[i-1];
    }

    STATE_VECTOR<TV>& operator=(const VECTOR_ND<T>& v)
    {
        assert(v.Size()==length);
        for(int i=1;i<=length;i++)
        {
            (*this)(i)=v(i);
        }
        return (*this);
    }

//#####################################################################
};

template<class TV> std::ostream& operator<<(std::ostream& output, const STATE_VECTOR<TV>& s)
{
    output << "[";
    for(int i=1;i<=s.Total_Size();i++)
    {
        output << s(i) << " ";
    }
    output << "]" << std::endl;
    return output;
}
}
#endif
