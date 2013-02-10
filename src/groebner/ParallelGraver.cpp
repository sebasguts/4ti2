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

#include "groebner/ParallelGraver.h"
#include "groebner/Globals.h"
#include "groebner/HermiteAlgorithm.h"
#include "groebner/CircuitMatrixAlgorithm.h"
#include "groebner/LatticeBasis.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/AugmentedVectorArray.h"

#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/Norms.hpp"

#include <algorithm>
#include <iostream>
#include <future>
#include <cstdlib>

using namespace _4ti2_;


ParallelGraver::ParallelGraver()
{
}

ParallelGraver::~ParallelGraver()
{
}

void print_vector (const std::vector<int>& v){
    for (auto it = v.begin(); it != v.end(); it++)
	*out << *it << " ";
    *out << "\n";
}


/** 
 * Compute Graver basis
 * 
 * @param feasible Feasible Problem
 * @param basis Vector array of the correct width to store result
 */
void
ParallelGraver::compute(
                Feasible& feasible,
                VectorArray& basis)
{
    // Some useful quantities
    int num_variables = feasible.get_dimension();

    *out << "Computing Graver basis (parallel) ...\n";

    // Step 1: Find a lattice basis
    lattice_basis (feasible.get_matrix(), basis);
    *out << "Lattice basis computed\n";
    *out << basis << std::endl;
    int m_rank = basis.get_number();

    *out << "Rank of problem: " << m_rank << "\n";

    // Step 2: Permute columns such that a full rank matrix is in the
    // first r columns.  P stores the permutation which undoes the
    // shuffling in the end.
    Permutation P = permute_full_rank_to_left (basis);

    // Project and compute a Graver basis with zsolve.
    VectorArray *projected = new VectorArray (basis.get_number(), m_rank);
    VectorArray::project(basis, 0, m_rank-1, *projected);

    *out << "Projected Lattice\n";
    *out << *projected << std::endl;

    // Call old zsolve Graver code to determine minimal elements.  
    _4ti2_state *m_state = new _4ti2_zsolve_::GraverAPI <IntegerType> ();
    // _4ti2_state_set_options(m_state, argv, argc);
    _4ti2_matrix *m_matrix = m_state->create_matrix(projected->get_number(),
						    projected->get_size(),
						    "lat");
    // Fill matrix using the GMP API
    for (int i = 0; i < projected->get_number(); i++) {
  	for (int j = 0; j < projected->get_size(); j++) {
#ifdef _4ti2_GMP_
  	    m_matrix->set_entry_mpz_class (i,j, (*projected)[i][j]);
#elif defined(_4ti2_INT64_)
	    m_matrix->set_entry_int64_t (i,j, (*projected)[i][j]);
#elif defined(_4ti2_INT32_)
	    m_matrix->set_entry_int32_t (i,j, (*projected)[i][j]);
#endif
	}
    }
    m_state->compute();
    _4ti2_matrix *m_result = m_state->get_matrix("zhom");

    // Clear projected and insert vectors from the result (which is a
    // Graver basis of the projected lattice)

    VectorArray *projected_g_basis = new VectorArray(m_result->get_num_rows(), m_rank);
    projected_g_basis->renumber(m_result->get_num_rows());
    for (int i = 0; i < m_result->get_num_rows(); i++)
	for (int j = 0; j < m_result->get_num_cols(); j++){
#ifdef _4ti2_GMP
	    m_result->get_entry_mpz_class (i,j, (*projected_g_basis)[i][j]);
#elif defined(_4ti2_INT64_)
	    m_result->get_entry_int64_t (i,j, (*projected_g_basis)[i][j]);
#elif defined(_4ti2_INT32_)
	    m_result->get_entry_int32_t (i,j, (*projected_g_basis)[i][j]);
#endif
	}
    delete m_state; // Clean up _4ti2_state
    
    *out << "Graver basis of the projected lattice: \n";
    *out << *projected_g_basis;
    *out << "\n";

    // TODO: Move this lifting code inside the loop.  In the end we
    // only want to lift one step at a time.

    // need to transpose projected since graver moves are a combination of rows!
    VectorArray *projected_transposed = new VectorArray (projected->get_size(), projected->get_size());
    VectorArray::transpose (*projected, *projected_transposed);
    VectorArray *lifted_graver = lift_with_basis((*projected_g_basis), *projected_transposed, basis);
    delete projected_transposed;
    
    // Add negatives
    int stop = lifted_graver->get_number();
    for (int i = 0; i < stop; i++)
	lifted_graver->insert( - (*lifted_graver)[i] );

    *out << "Graver basis property on rank many components:\n";
    *out << *lifted_graver << std::endl;
    *out << "\n";

    // The big lifting loop. varindex is the actual index of the next
    // variable.  For instance if we have a kx3 matrix of already
    // lifted stuff, then varindex would start out as "3" and all
    // norms to be computed in the following are over 0,1,2.
    for (int varindex = lifted_graver->get_size()-1; varindex < num_variables; varindex++ ) {
	*out << "Now lifting variable number " << varindex+1 << "\n";

	// Create auxilliary data
	auto current = new AugmentedVectorArray (std::move(*lifted_graver));
	current->createNormBST (varindex);

	IntegerType max_norm = current->maximum_norm ();
	*out << "Maximum norm on first components:" << max_norm << "\n";
	// Debug: ask various things about the normBST, spit it out, etc.

	typedef _4ti2_zsolve_::NormPair< IntegerType > NP;
	typedef std::vector< NP > NPvector;
	NPvector jobs;
	if (current->get_tree()->cbegin()->first == 1)
	    jobs.push_back( NP(1,1) );
	IntegerType completed_norm = 1;
	while (completed_norm < 2*max_norm) {
	    std::vector < std::future < VectorArray* > > m_futures;
	    for (auto it = jobs.begin(); it != jobs.end(); it++) {
		if ((*it).sum <= completed_norm+1) {
		    // do this job ...
		    // Check if there are vectors with this norm:
		    if ( current->get_tree()->find(it->first) != current->get_tree()->end() &&
			 current->get_tree()->find(it->second) != current->get_tree()->end()) {
			std::future < VectorArray* > fut = std::async(
			    std::launch::async, // Compiler can decide launch order
			    &ParallelGraver::graverJob,
			    this, // Since GraverJob is a non-static member, it needs 'this' as an argument
			    *(current->get_tree()->at(it->first)),
			    *(current->get_tree()->at(it->second)),
			    *(current->get_vectors()),
			    varindex
			    );
			// For debug purposes, wait for finish?
			// fut.wait();
			m_futures.push_back (std::move(fut));
		    }
		}
	    } // for over jobs
	    // Synchronize:
	    *out << "Waiting for sync ... ";
	    for (auto it = m_futures.begin(); it != m_futures.end(); ++it)
		it->wait();
	    *out << "Done.\n";
	    // Retrieve results
	    for (auto it = m_futures.begin(); it != m_futures.end(); ++it) {
		VectorArray *res = it->get();
		if (res->get_number() == 0) // nothing new
		    continue;
		IntegerType norm = (*res)[0].norm(varindex);  // TODO: Norm need not be computed, is clear from pair (r,s)
		*out << "Current norm: " << norm << "\n";
		// check if maximum increased
		if (max_norm < norm) {
		    *out << "New norm bound : " << norm << "\n";
		    max_norm = norm;
		}
		// Store new Graver elements
		current->get_vectors()->insert (*res);
		// Update normBST
		if (current->get_tree()->find(norm) == current->get_tree()->end()) {
		    // not found
		    current->get_tree()->insert( std::pair <IntegerType, VectorArray* >
						 (norm, res));
		    // Ownership of the VectorArray transferred to current!
		}
		else {
		    // found norm, append vectors
		    /// @TODO Move semantics
		    for (Index vc = 0; vc < res->get_number(); vc++)
			current->get_tree()->at(norm)->insert( (*res)[vc] );
		    delete res;
		}
	    }
	    // Clean up futures:
	    m_futures.clear();
	    
	    // Now all jobs with norm sum <= completed_norm+1 are done.
	    // Schedule new jobs:
	    completed_norm++;
	    // Add new jobs for each norm pair (i, completed_norm), i=
	    // 1..completed_norm such that there are moves in the
	    // respective degrees.
	    if (current->get_tree()->find(completed_norm) != current->get_tree()->end())
		for (IntegerType i = 1; i <= completed_norm; i++)
		    if (current->get_tree()->find(i) != current->get_tree()->end())
			jobs.push_back ( NP (i, completed_norm));
	    // Keep jobs sorted according to total norm
	    std::sort(jobs.begin(), jobs.end());
	}; // while (completed_norm < 2*max_norm)
	// Before stepping to the next lift, rebuild the NormBST
	std::cout << "Done with variable: " << varindex << " Updating norm map." << std::endl;
	if (varindex < num_variables )
	    current->createNormBST(varindex+1);
    } // big for loop over varindex


	
	
    } // End of the big variable lifting loop
        
    // Undo the permutation so that the users coordinates are
    // restored.
    basis.permute(P);
} // compute()

/** 
 * \brief Lift a Vector according to given lift on bases.
 * 
 * @param v Vector to be lifted
 * @param basis VectorArray consisting of the projections of lifted_basis
 * @param lifted_basis VectorArray of long basis vectors, aligned with basis
 * 
 * @return A new Vector that is the lift
 */
Vector* 
ParallelGraver::lift_with_basis( 
    const Vector& v,
    const VectorArray& basis,
    const VectorArray& lifted_basis)
{
    // use zsolve to solve (basis * x = v) (qsolve can only handle homogeneous systems atm)
    _4ti2_zsolve_::ZSolveAPI<IntegerType> *zsolve_api = new _4ti2_zsolve_::ZSolveAPI<IntegerType>;
    // How to pass options?
//    char *opt = "zsolve --quiet";
//    zsolve_api->set_options(2, &opt);
    std::stringstream connect_basis;
    std::stringstream connect_rhs;
    connect_basis << basis.get_number() << " " << basis.get_size() << " ";
    connect_basis << basis;
    connect_rhs << "1 " << v.get_size() << " ";
    connect_rhs << v;
    zsolve_api->create_matrix(connect_basis, "mat");
    zsolve_api->create_matrix(connect_rhs, "rhs");
    zsolve_api->compute();
    _4ti2_matrix *result = zsolve_api->get_matrix("zinhom"); // This is the desired coefficient vector
    Vector *coefficient_vector = new Vector (lifted_basis.get_number());
    for (int i =0; i < lifted_basis.get_number(); i++) {
#ifdef _4ti2_GMP
	    result->get_entry_mpz_class (0,i, (*coefficient_vector)[i]);
#elif defined(_4ti2_INT64_)
	    result->get_entry_int64_t (0,i, (*coefficient_vector)[i]);
#elif defined(_4ti2_INT32_)
	    result->get_entry_int32_t (0,i, (*coefficient_vector)[i]);
#elif 1
	    // This branch is only to remove -Wunused-variable warning
	    // that appears without it.  It should never be used
	    i = result->get_num_cols();
#endif
    }
    delete zsolve_api;
    
    Vector *result_vector = new Vector (lifted_basis.get_size(), 0);
    for (int i = 0; i < lifted_basis.get_number(); i++)
	result_vector->add( (lifted_basis[i]) * ((*coefficient_vector)[i]));
    delete coefficient_vector;
    return result_vector;
}

VectorArray*
ParallelGraver::lift_with_basis( 
    const VectorArray& va,
    const VectorArray& basis,
    const VectorArray& lifted_basis)
{
    VectorArray *result = new VectorArray (0, lifted_basis.get_size());
    for (int i = 0; i < va.get_number(); i++)
	// @TODO: std::move!
	result->insert( (*lift_with_basis (va[i], basis, lifted_basis)) );
    return result;
}


/** 
 * Permute a full rank matrix to the left
 * 
 * @param va VectorArray to be changed
 * 
 * @return a permutation that undoes the permutation done in the function.
 */
Permutation
ParallelGraver::permute_full_rank_to_left (VectorArray& va){
    // Algorithm: Compute Hermite form and identify columns where
    // pivot is in next row.
    VectorArray *temp = new VectorArray (va);
    hermite(*temp);
    std::vector<int> piv;
    std::vector<int> npiv;
    for (int i = 0; i < temp->get_size(); i++) {
	if (piv.size() == (uint) va.get_number())
	    // Have enough vectors already
	    npiv.push_back(i);
	else if ((*temp)[i-npiv.size()][i] != 0)
	    // Rank extending column found!
	    piv.push_back(i);
	else 
	    npiv.push_back(i);
    }
    delete temp;
    Permutation p;
    for (auto it = piv.begin(); it != piv.end(); it++)
	p.push_back(*it);
    for (auto it = npiv.begin(); it != npiv.end(); it++)
	p.push_back(*it);
    va.permute(p);
    // Now we invert the permutation p such that we can return the
    // inverse permutation
    Permutation pinv;
    for (int i = 0; i < va.get_size(); i++){
	// Where is i in p?
	int index_of_i = 0;
	while (p[index_of_i] != i) 
	    index_of_i++;
	pinv.push_back(index_of_i);
    }
    return pinv;
}


bool
is_reducible (const Vector& v) {
    return true;
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
 * 
 */
VectorArray*
ParallelGraver::graverJob (const VectorArray& Gr,
			   const VectorArray& Gs,
			   const VectorArray& current_gens,
			   const Index& maxvar) 
{
    VectorArray *result = new VectorArray (0,current_gens.get_size());
    Vector *v;
    for (int i = 0; i < Gr.get_number(); i++ ) {
	for (int j = 0; j < Gs.get_number(); j++ ) {
	    v = new Vector (Gr[i] + Gs[j]);
	    if (!is_reducible (*v))
		result->insert (std::move(*v));
	    else 
		delete v;
	}
    }
    return result;
}
