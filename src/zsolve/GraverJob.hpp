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

#ifndef _4ti2_zsolve__GraverJob_
#define _4ti2_zsolve__GraverJob_

#include <iostream>
#include <cassert>

#include "Lattice.hpp"

namespace _4ti2_zsolve_
{

/**
 * \brief Graver job computation
 * 
 * Input:
 *   * lattice basis L for \Lattice\subseteq\Z^n (typically, a coordinate projection of ker(A))
 *   * sets of vectors G_r,G_s\in\Lattice with 1-norms r,s on first d components
 *   * the auxiliary information for each vector should be given
 *   from master to slave!
 *   * d
 *   * n = number of variables/components (Typically, we have d=n-1.)
 * 
 * Output:
 *   * all sums v+w \in G_r+G_s that are minimal elements in \Lattice\setminus\{0\}
 * 
 * Algorithm:
 *   * irr={}
 *   * for all pairs (v,w) \in G_r+G_s {
 *     * if v+w is minimal in \Lattice\setminus\{0\} {
 *       * add auxiliary data to v+w like support, norm, ...
 *       * set irr=irr\cup{v+w}
 *     }
 *   }
 *   * return irr (including auxiliary data for each vector)
 * 
 * 
 * Remarks:
 * 
 * * Clearly, if r=s, only half the number of pairs have to be
 * checked.
 * * If v and w have the same sign pattern in the last n-d
 * components, the quick test applies and v+w is reducibly by both
 * v and w and can thus be discarded immediately.
 * * Consequently, if d=n-1, it is worth pre-sorting G_r and G_s
 * according to the last component being >0, =0, <0.
 * * For the minimality test, we could use the lattice basis and
 * look for a nonzero vector u<>v+w in \Lattice with the same
 * sign-pattern as v+w and with u<=v+w in this orthant.
 * * We may also give the slave G_1, G_2, ...up to some norm for a
 * quick check.
 * 
 * 
 */

// template <typename T>
// VectorArray<T>*
// graverJob () {
//     return new VectorArray<T>;
// }

template <typename T>
bool
is_reducible (T* v) {
    return true;
}

template <typename T> 
VectorArray<T> *
graverJob (const VectorArray<T>& Gr,
	   const VectorArray<T>& Gs,
	   const Lattice<T>& current_gens,
	   const size_t& maxvar) 
{
    size_t num_vars = current_gens.num_variables();
    VectorArray<T> *result = new VectorArray<T> (num_vars);

    T* v = create_vector<T> (num_vars);
    for (size_t i = 0; i < Gr.num_vectors(); i++ ) {
	for (size_t j = 0; j < Gs.num_vectors(); j++ ) {
	    add_vectors (Gr[i], Gs[j], v, num_vars);
	    if (!is_reducible (v))
		result->append_vector ( copy_vector<T> (v, num_vars));
	}
    }
    return result;
}

} // namespace _4ti2_zsolve_

#endif
