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

#include <vector>
#include <iostream>

#include "groebner/GraverVector.h"
#include "groebner/BitSet.h"

using namespace _4ti2_;

void 
GraverVector::fill_supports_and_norm() {
    norm = 0;
    // norm is on the first n-1 entries only
    for (int i = 0; i < v->get_size()-1; i++){
	if ((*v)[i] > 0) {
	    pos.set(i);
	    norm += (*v)[i];
	}
	else if ((*v)[i] < 0) {
	    neg.set(i);
	    norm -= (*v)[i];
	}
    }
    if ((*v)[v->get_size()-1] > 0)
	pos.set(v->get_size()-1);
    if ((*v)[v->get_size()-1] < 0)
	neg.set(v->get_size()-1);
}


GraverVector::GraverVector (Vector* _v) : pos (_v->get_size(), 0), neg (_v->get_size(), 0), v(_v)
{
    fill_supports_and_norm();
}

GraverVector::GraverVector (const GraverVector& _v) : pos (_v.pos), neg (_v.neg), norm (_v.norm), v (_v.v)
{ 
    ;
}
									 
GraverVector::GraverVector (GraverVector&& _v) : pos (_v.pos), neg (_v.neg), norm (_v.norm), v (std::move (_v.v))
{
    ;
}

void
GraverVector::swap (GraverVector & other) {
    std::swap(v, other.v);
    std::swap(norm, other.norm);
    std::swap(pos, other.pos);
    std::swap(neg, other.neg);
}

// GraverVector& 
// GraverVector::operator=(GraverVector other)
// {
//     swap (other); 
//     return *this;
// }

GraverVector& 
GraverVector::operator=(const GraverVector& other)
{
    v = other.v;
    pos = other.pos;
    neg = other.neg;
    norm = other.norm;
    return *this;
}
    
GraverVector& 
GraverVector::operator=(GraverVector&& other) {
    v = std::move (other.v);
    pos = other.pos;
    neg = other.neg;
    norm = other.norm;
    return *this;
}

GraverVector::~GraverVector() {

}
