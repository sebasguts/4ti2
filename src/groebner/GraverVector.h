/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2013 4ti2 team.
Main author(s): Thomas Kahle.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#ifndef _4ti2_groebner__GraverVector_
#define _4ti2_groebner__GraverVector_

#include "groebner/BitSet.h"
#include "groebner/Vector.h"
#include "groebner/Size.h"

namespace _4ti2_
{

// A graver Vector is an overlay to an actual vector.  It stores
// supplementary data together with a pointer to the actual vector.
class GraverVector
{
public:
    GraverVector (Vector* _v);
    GraverVector (const GraverVector& _v);
    GraverVector (GraverVector&& _v);
    // GraverVector& operator=(GraverVector);
    GraverVector& operator=(const GraverVector&);
    GraverVector& operator=(GraverVector&&);
    ~GraverVector();

    void swap (GraverVector& other);
    void fill_supports_and_norm();

    IntegerType last_entry () const { return (*v)[v->get_size()-1];}; // Todo: Optimize

    bool is_sign_consistent (const GraverVector& other) const;
    bool is_sign_support_below (const GraverVector& other) const;
    bool is_applicable (const GraverVector& other) const;

    bool is_below (const GraverVector& other) const;

    BitSet pos;
    BitSet neg;
    IntegerType norm; ///< Norm on first n-1 components!
    Vector* v;
};

inline
bool
GraverVector::is_sign_consistent (const GraverVector& other) const {
    return ( BitSet::set_disjoint (neg, other.pos) &&
	     BitSet::set_disjoint (pos, other.neg) );
}

inline
bool
GraverVector::is_below (const GraverVector& other) const {
    assert (v->get_size() == other.v->get_size());
    for (int i = 0; i < v->get_size(); i++){
	if ((*v)[i] > 0 && (*v)[i] > (*other.v)[i])
	    return false;
	if ((*v)[i] < 0 && (*v)[i] < (*other.v)[i])
	    return false;
    }
    return true;
}
inline
bool
GraverVector::is_sign_support_below (const GraverVector& other) const {
    return ( BitSet::set_subset (pos, other.pos) &&
	     BitSet::set_subset (neg, other.neg) );
}

inline
bool
GraverVector::is_applicable (const GraverVector& other) const {
    // sign consistent everywhere except in the last position
    const Size size = v->get_size();
    BitSet n (neg);
    BitSet p (pos);
    n.flip (size-1);
    p.flip (size-1);
    return ( BitSet::set_disjoint (n, other.pos) &&
	     BitSet::set_disjoint (p, other.neg) );
}

} // namespace _4ti2_

#endif
