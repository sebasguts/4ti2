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

#ifndef _4ti2_groebner__ReductionTree_
#define _4ti2_groebner__ReductionTree_

#include <map>

#include "groebner/GraverVectors.h"
#include "groebner/GraverTypes.h"

namespace _4ti2_
{

class ReductionTree;
typedef std::map<IntegerType, ReductionTree> TreeMap;

class ReductionTree
{
public:
    ReductionTree ();
    ~ReductionTree ();

    void clear (); // Reset the tree;
    void insert (const Vector& v); ///< Insert a new vector into the tree structure
    bool isReducible (const Vector& v) const;  ///< Returns a wether given vector is reducible modulo the stored vectors

private:
    void insert_with_offset (const Vector& v, const int offset);
    bool isReducible_with_offset (const Vector& v, const int offset) const; 
    TreeMap m_branches;
    // GraverVector *m_vector; // Probably we never need to produce the reducing vector
};

} // namespace _4ti2_

#endif
