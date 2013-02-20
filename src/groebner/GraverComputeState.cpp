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

#include "groebner/GraverComputeState.h"
#include "groebner/GraverVectors.h"
#include "groebner/GraverVectorsNaive.h"
#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"

#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/GraverAPI.hpp"

using namespace _4ti2_;

typedef GraverVectorsNaive GraVec;

GraverComputeState::GraverComputeState (const VectorArray& lb){
    m_graverVectors = new GraVec (lb);
    m_latticeBasis = new VectorArray (lb);
    m_rank = lb.get_number(); // Trust that we got an actual basis!

    for (int i = 0; i < lb.get_size(); i++){
	m_projected_lattice_bases.push_back (new VectorArray(lb.get_number(), i));
	VectorArray::project(lb, 0, i, *m_projected_lattice_bases.back());
	// std::cout << "Projected basis: " << *m_projected_lattice_bases.back();
    }
    // Insert the full basis too for lifting convenience
    m_projected_lattice_bases.push_back(new VectorArray (lb));
}

GraverComputeState::~GraverComputeState() {
    delete m_graverVectors;
    delete m_latticeBasis;
    for (auto it = m_projected_lattice_bases.begin(); it != m_projected_lattice_bases.end(); ++it)
	delete *it;
}

VectorArray
GraverComputeState::get_vectors () {
    VectorArray v(0,m_graverVectors->get_size());
    for (auto it = m_graverVectors->begin(); it != m_graverVectors->end(); ++it){
	v.insert(*it);
    }
    return v;
}


VectorArray 
GraverComputeState::get_vectors_without_negatives () {
    VectorArray v = get_vectors();
    VectorArray result (0,v.get_size());
    for (int i = 0; i < v.get_number(); i++) {
	bool has_negative = false;
	for (int j = i+1; j < v.get_number(); j++) {
	    if ( is_below ( -(v[i]), v[j] ) ) {
		has_negative = true;
		break;
	    }
	}
	if (!has_negative)
	    (v[i] <= -v[i]) ? result.insert( std::move(-v[i]) ) : result.insert(std::move(v[i]));
    }
    return result;
}

void 
GraverComputeState::projectToRank () {
    delete m_graverVectors;
    m_graverVectors = new GraVec (*m_projected_lattice_bases[m_rank-1]);
}

inline void
GraverComputeState::MakeGraverBasisWithZSolve() 
{
    VectorArray v( get_vectors() );
    MakeGraverBasisWithZSolve (v);
    delete m_graverVectors;
    m_graverVectors = new GraVec ( std::move(v) );
};

void 
GraverComputeState::liftOneStep () {
    m_graverVectors->lift (*m_projected_lattice_bases [m_graverVectors->get_size()+1]);
}

void 
GraverComputeState::liftGraverProperty () {
    // Cheating
    MakeGraverBasisWithZSolve();    
}

#if 0
IntegerType 
GraverComputeState::maximum_norm () const {
    assert(m_normTree != NULL);
    // std::cout << m_normTree->rbegin()->first << "\n";
    return m_normTree->rbegin()->first;
}
#endif

/** 
 * Extend a VectorArray to be a GraverBasis (without negatives)
 * 
 * @param basis lattice basis to be extended to a Graver basis
 */
void
GraverComputeState::MakeGraverBasisWithZSolve (VectorArray& basis){
    // Call old zsolve Graver code to determine minimal elements.  
    _4ti2_state *m_state = new _4ti2_zsolve_::GraverAPI <IntegerType> ();
    // _4ti2_state_set_options(m_state, argv, argc);
    _4ti2_matrix *m_matrix = m_state->create_matrix(basis.get_number(),
						    basis.get_size(),
						    "lat");
    // Fill matrix using the GMP API
    for (int i = 0; i < basis.get_number(); i++) {
  	for (int j = 0; j < basis.get_size(); j++) {
#ifdef _4ti2_GMP_
  	    m_matrix->set_entry_mpz_class (i,j, basis[i][j]);
#elif defined(_4ti2_INT64_)
	    m_matrix->set_entry_int64_t (i,j, basis[i][j]);
#elif defined(_4ti2_INT32_)
	    m_matrix->set_entry_int32_t (i,j, basis[i][j]);
#endif
	}
    }
    m_state->compute();
    _4ti2_matrix *m_result = m_state->get_matrix("zhom");
    /// Potentially do better here if profiling shows that this matters
    /// Problem: We write information to basis that is already there.
    /// Unsure if we can savely skip the first rows because the called
    /// code may have permuted the rows?
    basis.renumber(m_result->get_num_rows());
    for (int i = 0; i < m_result->get_num_rows(); i++)
	for (int j = 0; j < m_result->get_num_cols(); j++){
#ifdef _4ti2_GMP
	    m_result->get_entry_mpz_class (i,j, basis[i][j]);
#elif defined(_4ti2_INT64_)
	    m_result->get_entry_int64_t (i,j, basis[i][j]);
#elif defined(_4ti2_INT32_)
	    m_result->get_entry_int32_t (i,j, basis[i][j]);
#endif
	}
    delete m_state; // Clean up _4ti2_state

    // Add negatives because the zsolve Graver basis does not include
    // those
    int stop = basis.get_number();
    for (int i = 0; i < stop; i++){
	basis.insert(-basis[i]);
    }
}

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
 */


VectorArray*
GraverComputeState::graverJob (const VectorArray& Gr,
			   const VectorArray& Gs,
			   const VectorArray& current_gens,
			   const Index& maxvar)
// TODO: What is maxvar good for?
{
    // Step 1: Compute all sums f + g where f and g are sign
    // consistent with f \in Gr, g \in Gs.
    VectorArray *result = new VectorArray (0,Gr.get_size());
    for (int i = 0; i < Gr.get_number(); i++ ) {
	// if Gr and Gs are the same, then we don't need all pairs:
	int j = 0;
	if (&Gr == &Gs) 
	    j = i+1;
	for (; j < Gs.get_number(); j++ ) {
	    // Check for sign consistency on the first n-1 components
	    if (
		(sign_consistent(Gr[i], Gs[j], Gr.get_size()-1))
		&&
		// not sign consistent on last variable
		( 
		    (Gr[i][maxvar] <0 && Gs[j][maxvar] > 0)
		    || 
		    (Gr[i][maxvar] > 0 && Gs[j][maxvar] < 0)
		    )
		) {
		Vector sum = Gr[i] + Gs[j];
		std::cout << "I decided to check the sum " << Gr[i] << " + " << Gs[j] << "=" << sum << "\n" ;
		if (GraverComputeState::is_reducible(sum, current_gens)) {
		    std::cout << "oops, was reducible... \n";
		}
		else {
		    std::cout << "great, not reducible ! \n";
		    result->insert(Vector(sum));
		}
		// result->insert (new Vector (Gr[i] + Gs[j]));
	    }
	}
    }
    return result;
}

#if 0
void
GraverComputeState::createNormBST (Index stop) {
    std::cout << "Creating norm tree for indices 0 .. " << stop-1 << std::endl;
    delete m_normTree;
    m_stopindex = stop;
    m_normTree = new NormBST ();
    /// @TODO: Parallelize this computation by splitting the todolist
    for (int i = 0; i < m_array->get_number(); i++){
	IntegerType current_norm = (*m_array)[i].norm(stop);
	VectorArray *current_vectors = NULL;
	auto it = m_normTree->find(current_norm);
	if (it != m_normTree->end()) {
	    current_vectors = it->second;
	    current_vectors->insert( (*m_array)[i] );
	} 
	else {
	    current_vectors = new VectorArray (0,m_array->get_size());
	    current_vectors->insert( (*m_array)[i] );
	    m_normTree->insert( std::pair <IntegerType, VectorArray* > (current_norm, current_vectors) );
	};
    }
    std::cout << "Norm Tree created, minimum norm: " << m_normTree->begin()->first;
    std::cout << ", maximum norm : " << m_normTree->rbegin()->first << "\n";
}
#endif 
