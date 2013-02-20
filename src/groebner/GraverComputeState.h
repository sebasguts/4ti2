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

#ifndef _4ti2_groebner__GraverComputeState_
#define _4ti2_groebner__GraverComputeState_

#include <vector>
#include <map>

#include "groebner/GraverVectors.h"
#include "groebner/VectorArray.h"

namespace _4ti2_
{

class GraverComputeState
{
public:
    GraverComputeState(const VectorArray& lb); ///< copy lattice basis
    ~GraverComputeState();

    VectorArray get_vectors (); ///< Current state of computation
    VectorArray get_vectors_without_negatives (); ///< Get vectors without duplicating negatives

    void projectToRank (); ///< Project down to rank of the problem
    void MakeGraverBasisWithZSolve ();

    void liftOneStep (); ///< Lift one dimension
    void liftGraverProperty (); ///< Make the recently lifted Vectors a GraverBasis

private:
    static void MakeGraverBasisWithZSolve (VectorArray& basis);

    static bool is_below (const Vector& v1, const Vector& v2);
    static bool is_reducible (const Vector& v, const VectorArray& va);
    static bool sign_consistent (const Vector& v1, const Vector& v2, Index stop);

    static VectorArray* graverJob (const VectorArray& Gr,
				   const VectorArray& Gs,
				   const VectorArray& current_gens,
				   const Index& maxvar);

// Data:
private:
    VectorArray *m_latticeBasis;
    GraverVectors *m_graverVectors;
    std::vector < VectorArray* > m_projected_lattice_bases;

    // NormBST *m_normTree;
    // Index m_stopindex; //?
    Size m_rank;
};

inline bool
GraverComputeState::is_below (const Vector& v1, const Vector& v2) {
    assert (v1.get_size() == v2.get_size());
    for (int i = 0; i< v1.get_size(); i++){
	if (v1[i] > 0 && v1[i] > v2[i])
	    return false;
	if (v1[i] < 0 && v1[i] < v2[i])
	    return false;
    }
    return true;
}

inline bool 
GraverComputeState::sign_consistent (const Vector& v1, const Vector& v2, Index stop)
{
    for (Index i = 0; i < stop; i++) {
	if (v1[i] <0 && v2[i] >0 )
	    return false;
	if (v1[i] > 0 && v2[i] < 0 )
	    return false;
    }
    return true;
}

inline bool
GraverComputeState::is_reducible (const Vector& v, const VectorArray& va) {
    assert (v.get_size() == va.get_size());
    for (int i=0; i<va.get_number(); i++)
	if (GraverComputeState::is_below (va[i], v))
	    return true;
    return false;
}


} // namespace _4ti2_

#endif
