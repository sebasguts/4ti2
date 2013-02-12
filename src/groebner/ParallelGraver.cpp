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
#include "groebner/GraverComputeState.h"

#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/GraverAPI.hpp"
#include "zsolve/Norms.hpp"

#include <algorithm>
#include <iostream>
#include <future>
#include <cstdlib>

#include <system_error>

using namespace _4ti2_;


ParallelGraver::ParallelGraver()
{
}

ParallelGraver::~ParallelGraver()
{
}

void print_mark() {
    static int i = 0;
    std::cout << "Mark no. "<< i++ << "\n";
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

    GraverComputeState current (basis);  // give lattice basis for initialization

    // Compute a Graver basis of the smallest projection with zsolve.
    current.get_vectors() = current.get_projected_basis(0);
    MakeGraverBasisWithZSolve (current.get_vectors());
    
    *out << "Graver basis of the projected lattice: \n";
    *out << current.get_vectors();
    *out << "\n";

    // Add negatives because the zsolve Graver basis does not include
    // those
    int stop = current.get_vectors().get_number();
    for (int i = 0; i < stop; i++){
	current.get_vectors().insert(-current.get_vectors()[i]);
    }

    // Graver state is initialized with a negative containing Graver
    // basis now.  Time to start the lifting.

    // The big lifting loop. varindex is the actual index of the next
    // variable.  For instance if we have a kx3 matrix of already
    // lifted stuff, then varindex would start out as "3" and all
    // norms to be computed in the following are over 0,1,2.
    for (Index varindex = m_rank; varindex < num_variables; varindex++ ) {
	*out << "Now lifting variable number " << varindex+1 << "\n";

	current.liftToIndex (varindex-m_rank);
    
// 	*out << "Graver basis property on rank many components:\n";
// 	*out << current.get_vectors() << std::endl;
// 	*out << "\n";

	// Create auxilliary data
	current.createNormBST (varindex);

	IntegerType max_norm = current.maximum_norm ();
	*out << "Maximum norm on first components:" << max_norm << "\n";
	// Debug: ask various things about the normBST, spit it out, etc.

	typedef _4ti2_zsolve_::NormPair< IntegerType > NP;
	typedef std::vector< NP > NPvector;
	NPvector jobs;
	if (current.get_tree().cbegin()->first == 1)
	    jobs.push_back( NP(1,1) );
	IntegerType current_norm = 1;
	while (current_norm < 2*max_norm) {
	    std::cout << "Now doing norm: "<< current_norm << "\n";
	    // std::vector < VectorArray * > m_futures; // The future is now!
	    std::vector < std::future < VectorArray* > > m_futures;  // The future is later
	    // find lowest norm jobs:
	    NPvector lowest_norm_jobs;
	    std::copy_if (jobs.begin(),
			  jobs.end(),
			  std::back_inserter(lowest_norm_jobs),
			  [current_norm](const NP& np) -> bool 
			  { return (np.sum == current_norm);} );
	    auto new_end = std::remove_if (jobs.begin(),
					   jobs.end(),
					   [current_norm](const NP& np) -> bool 
					   { return (np.sum == current_norm);} );
 	    jobs.erase (new_end, jobs.end());
 	    std::cout << "Lowest norm jobs : ";
 	    for (auto it = lowest_norm_jobs.begin(); it != lowest_norm_jobs.end(); it++) {
 		std::cout << "(" << it->first << "," << it->second << "),";
 	    }
	    std::cout << "\n";
 	    std::cout << "Other jobs: ";
 	    for (auto it = jobs.begin(); it != jobs.end(); it++) {
 		std::cout << "(" << it->first << "," << it->second << "),";
 	    }
 	    std::cout << "\n";
	    if (lowest_norm_jobs.size() == 0) {
		std::cout << "No jobs of norm "<< current_norm << "\n";
	    }
	    // Do lowest norm jobs:
	    for (auto it = lowest_norm_jobs.begin(); it != lowest_norm_jobs.end(); it++) {
		std::cout << "Doing job :" << it->first << "," << it->second << "\n";
		// Check if there are vectors with this norm:
		if ( current.get_tree().find(it->first) != current.get_tree().end() &&
		     current.get_tree().find(it->second) != current.get_tree().end()) {
// 		    VectorArray *fut = graverJob (
// 			*(current.get_tree().at(it->first)),
// 			*(current.get_tree().at(it->second)),
// 			current.get_vectors(),
// 			varindex
// 			);
		    try {
			std::future < VectorArray* > fut = std::async(
			    std::launch::async, // This directive makes it launch a new thread for each job (not good if there are many!)
			    ParallelGraver::graverJob,
			    *(current.get_tree().at(it->first)),
			    *(current.get_tree().at(it->second)),
			    current.get_vectors(),
			    varindex
			    );
			m_futures.push_back (std::move(fut));
		    }
		    catch (std::system_error e) {
			std::cout << ":(\n";
			std::exit (1);
		    }
		    
		    // For debug purposes, wait for finish?
		    // fut.wait();
		} // if (current. ...
	    } // for over lowest_degree_jobs
	    // Synchronize:
 	    *out << "Waiting for sync ... ";
 	    for (auto it = m_futures.begin(); it != m_futures.end(); ++it)
 		it->wait();
 	    *out << "Done.\n";
// 	    // Retrieve results
	    for (auto it = m_futures.begin(); it != m_futures.end(); ++it) {
		VectorArray *res = it->get();
//		VectorArray *res = *it;
		if (res->get_number() == 0) // nothing new
		    continue;
//		IntegerType norm = (*res)[0].norm(varindex);  // TODO: Norm need not be computed, is clear from pair (r,s)
		*out << "Current norm: " << current_norm << "\n";
		// check if maximum increased
		if (max_norm < current_norm) {
		    *out << "New norm bound : " << current_norm << "\n";
		    max_norm = current_norm;
		}
		std::cout << "I'm going to add new vectors. So far I got " << current.get_vectors().get_number() << std::endl;
		// Store new Graver elements
		current.get_vectors().insert (*res);
		std::cout << "and now there are: " << current.get_vectors().get_number() << std::endl;
		// Update normBST
		if (current.get_tree().find(current_norm) == current.get_tree().end()) {
		    // not found
		    current.get_tree().insert( std::pair <IntegerType, VectorArray* >
						 (current_norm, res));
		    // Ownership of the VectorArray transferred to current!
		}
		else {
		    // found norm, append vectors
		    current.get_tree().at(current_norm)->insert( std::move(*res));
		}
	    }
	    // Clean up futures:
	    m_futures.clear();
	    
	    // Now all jobs with norm sum <= current_norm+1 are done.
	    // Schedule new jobs:
	    current_norm++;
	    // Add new jobs for each norm pair (i, current_norm), i=
	    // 1..current_norm such that there are moves in the
	    // respective degrees.
	    std::cout << "Current Jobs overview: \n";
	    for (auto it=jobs.begin(); it != jobs.end(); it++)
		std::cout <<"("<<it->first<<","<<it->second<<"), ";
	    std::cout << std::endl;
	    std::cout << "Current size of Graver basis: " << current.get_vectors().get_number() << "\n";
	    if (current.get_tree().find(current_norm) != current.get_tree().end())
		for (IntegerType i = 1; i <= current_norm; i++)
		    if (current.get_tree().find(i) != current.get_tree().end())
			jobs.push_back ( NP (i, current_norm));
	    // Keep jobs sorted according to total norm
	    std::sort(jobs.begin(), jobs.end());
	}; // while (current_norm < 2*max_norm)
	// Before stepping to the next lift, rebuild the NormBST
	std::cout << "Done with variable: " << varindex << std::endl;
	if (varindex+1 < num_variables ) {
	    std::cout << "Updating norm map for next variable." << std::endl;
	    current.createNormBST(varindex+1);
	}
	std::cout << "Ok buddy, I'm done with this variable, here's what I got\n";
	std::cout << current.get_vectors();
	// TODO: At this point "current_graver_basis" will be destructed because the destructor of the 
    } // End of the big variable lifting loop
    // Undo the permutation so that the users coordinates are
    // restored.
    current.get_vectors().permute(P);
    // Todo: Do the reduction in a smarter way:
    std::cout << "Removing negatives...";
    basis.clear();
    for (int i = 0; i < current.get_vectors().get_number(); i++) {
	bool has_negative = false;
	for (int j = i+1; j< current.get_vectors().get_number(); j++) {
	    if ( ParallelGraver::is_below ( -(current.get_vectors()[i]), current.get_vectors()[j] ) ) {
		has_negative = true;
		break;
	    }
	}
	if (!has_negative)
	    basis.insert(Vector(current.get_vectors()[i]));
    }
    std::cout << " done\n";
    basis.sort();
    
} // compute()

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
ParallelGraver::is_reducible (const Vector& v, const VectorArray& va) {
    assert (v.get_size() == va.get_size());
    for (int i=0; i<va.get_number(); i++)
	if (ParallelGraver::is_below (va[i], v))
	    return true;
    return false;
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
ParallelGraver::graverJob (const VectorArray& Gr,
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
//		std::cout << "I decided to check the sum " << Gr[i] << " + " << Gs[j] << "=" << sum << "\n" ;
		if (ParallelGraver::is_reducible(sum, current_gens)) {
//		    std::cout << "oops, was reducible... \n";
		}
		else {
//		    std::cout << "great, not reducible ! \n";
		    result->insert(Vector(sum));
		}
		// result->insert (new Vector (Gr[i] + Gs[j]));
	    }
	}
    }
    return result;
}

void 
ParallelGraver::MakeGraverBasisWithZSolve (VectorArray& basis){
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
}
