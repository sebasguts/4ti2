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

#ifndef _4ti2_groebner__GraverVectorsNaive_
#define _4ti2_groebner__GraverVectorsNaive_

#include "groebner/VectorArray.h"
#include "groebner/Vector.h"
#include "groebner/GraverVectors.h"

namespace _4ti2_
{

/** 
 * Naive implementation of Graver Vectors: Pretty plain list
 */
class GraverVectorsNaive : public GraverVectors
{
public:
    GraverVectorsNaive ();
    GraverVectorsNaive (const VectorArray& va); ///< Copy a VectorArray
    ~GraverVectorsNaive();

    // Iterators to access the vectors
    GraverVectorsIterator begin () {return m_data2.begin();};
    GraverVectorsIterator end () {return m_data2.end();};
    GraverVectorsIterator begin (IntegerType norm); ///< Access by norm on first n-1 components
    GraverVectorsIterator end (IntegerType norm);  ///< Access by norm on first n-1 components

    void addNegatives ();  ///<  Insert the negative for each vector (if not already present)
    void removeNegatives (bool lexicographic = true); ///< Keep only v or -v for each vector

    Size get_size(); ///< Dimension aka length aka size
    
    void insert (Vector v); ///< Insert a new Vector
    void insert (Vector&& v); ///< Insert a new Vector

    void lift (const VectorArray& lifted_basis); ///< Lift to next dimension

private:
    // Lifting helpers
    Vector lift_with_basis( const Vector& v,
			    const VectorArray& basis,
			    const VectorArray& lifted_basis);

private:
    VectorArray *m_data;
    std::vector<Vector> m_data2;
    Size size;
};

} // namespace _4ti2_

#endif
