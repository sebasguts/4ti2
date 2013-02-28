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
#include "groebner/GraverTypes.h"
#include "groebner/ReductionTree.h"

namespace _4ti2_
{

/** 
 * Naive implementation of Graver Vectors: Pretty plain list
 *
 * Idea: The vectors in the array are stored in a big vector array.
 * Overlayed over this vector array are various containers storing
 * pointers to elements of the vector array.  For instance to access
 * all elements of a given norm one would get a std::vector of
 * pointers to elements in the VectorArray.
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
    Size get_number() { return m_data->get_number();}; ///< Number of Vectors stored

    VectorArray& get_vectors();
    IntegerType get_max_norm(); // only on the first n-1 elements

    // Hmm, the signature is different... that's a bit bad style
    VecVecP& get_vectors(IntegerType n);

    bool has_vectors_with_norm (IntegerType n);

    bool is_reducible (const Vector& v) const;

    void insert (Vector&& v); ///< Insert a new Vector (always move)
    void insert (VectorArray&& va);

    void lift (const VectorArray& lifted_basis); ///< Lift to next dimension

private:
    // Lifting helpers
    Vector lift_with_basis( const Vector& v,
			    const VectorArray& basis,
			    const VectorArray& lifted_basis);

    void createNormOL (); ///< Create the norm overlay from scratch
    void createRedTree (); ///< Create the reduction tree from scratch

    void check_sizes();

private:
    VectorArray *m_data;
    
    NormOL m_normOL; ///< Norm overlay
    ReductionTree m_redTree;

    std::vector<Vector> m_data2; // Dummy to satisfy the iterators
//    Size size;
};


inline
VectorArray&
GraverVectorsNaive::get_vectors() {
    return *m_data;
}

inline
VecVecP& 
GraverVectorsNaive::get_vectors(IntegerType n){
    assert (m_normOL.find(n) != m_normOL.end());
    return m_normOL.at(n);
}

inline
IntegerType 
GraverVectorsNaive::get_max_norm() {
    return m_normOL.rbegin()->first;
}

inline
bool
GraverVectorsNaive::has_vectors_with_norm (IntegerType n) {
    return !(m_normOL.find(n) == m_normOL.end());
}

inline
bool 
GraverVectorsNaive::is_reducible (const Vector& v) const {
    return m_redTree.isReducible (v);
}



} // namespace _4ti2_

#endif
