
#ifndef BOOST_MPL_SET_AUX_BEGIN_END_IMPL_HPP_INCLUDED
#define BOOST_MPL_SET_AUX_BEGIN_END_IMPL_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2003-2004
// Copyright David Abrahams 2003-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Source: /physbam_repository/External_Libraries/Archives/boost/boost/mpl/set/aux_/begin_end_impl.hpp,v $
// $Date: 2007/02/12 18:25:38 $
// $Revision: 1.1 $

#include <boost/mpl/begin_end_fwd.hpp>
#include <boost/mpl/set/aux_/iterator.hpp>

namespace boost { namespace mpl {

template<>
struct begin_impl< aux::set_tag >
{
    template< typename Set > struct apply
    {
        typedef s_iter<Set,Set> type;
    };
};

template<>
struct end_impl< aux::set_tag >
{
    template< typename Set > struct apply
    {
        typedef s_iter< Set,set0<> > type;
    };
};

}}

#endif // BOOST_MPL_SET_AUX_BEGIN_END_IMPL_HPP_INCLUDED
