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

#ifndef _4ti2_groebner__GraverVectors_
#define _4ti2_groebner__GraverVectors_

#include "groebner/Size.h"
#include "groebner/VectorArray.h"
#include "groebner/Vector.h"
#include "groebner/GraverVectors.h"
#include "groebner/GraverTypes.h"

namespace _4ti2_
{

/** 
 * Abstraction of a bunch of Graver Vectors.  Interface only.
 * 
 * @return 
 */
class GraverVectors
{
public:
    GraverVectors ();
    GraverVectors (const VectorArray& va); ///< Copy a VectorArray
    ~GraverVectors();

    virtual void addNegatives () = 0;  ///<  Insert the negative for each vector (if not already present)
    virtual void removeNegatives (bool lexicographic = true) = 0; ///< Keep only v or -v for each vector

    virtual Size get_size() = 0; ///< Dimension aka length aka size
    virtual Size get_number() = 0; ///< Number of stored vectors

    virtual VectorArray& get_vectors() = 0;
    virtual IntegerType get_max_norm() = 0; // only on the first n-1 elements

    virtual VecVecP& get_vectors(IntegerType n) = 0;

    // Access by norm and positive and negative support
    virtual VecVecP& get_vectors(IntegerType n, std::pair <LongDenseIndexSet, LongDenseIndexSet>) = 0;

    virtual bool has_vectors_with_norm (IntegerType n) = 0;

    virtual bool is_reducible (const Vector& v) const = 0;
    
//    virtual void insert (Vector v) = 0; ///< Insert a new Vector
    virtual void insert (Vector&& v) = 0; ///< Insert a new Vector
//    virtual void insert (VectorArray va) = 0;
    virtual void insert (VectorArray&& va) = 0;

    virtual void lift (const VectorArray& lifted_basis) = 0; ///< Lift to next dimension
};

} // namespace _4ti2_

#endif
