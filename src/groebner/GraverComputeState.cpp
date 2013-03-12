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

#include <algorithm>
#include <future>
#include <map>
#include <vector>
#include <iostream>

#include "groebner/BitSetStream.h"
#include "groebner/GraverComputeState.h"
#include "groebner/GraverTypes.h"
#include "groebner/GraverVectorsWithFilter.h"
#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"

#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/GraverAPI.hpp"

using namespace _4ti2_;

GraverComputeState::GraverComputeState (const VectorArray& lb){
    // the false call saves creation of aux data on the first call.
    m_graverVectors = new GraVec (lb, false);
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

inline bool
is_below (const Vector& v1, const Vector& v2) {
    assert (v1.get_size() == v2.get_size());
    for (int i = 0; i< v1.get_size(); i++){
	if (v1[i] > 0 && v1[i] > v2[i])
	    return false;
	if (v1[i] < 0 && v1[i] < v2[i])
	    return false;
    }
    return true;
}

VectorArray
GraverComputeState::get_vectors_without_negatives_destructive () {
    return m_graverVectors->get_vectors_without_negatives_destructive();
}

void
GraverComputeState::projectToRank () {
    std::cout << "In projectToRank\n";
    delete m_graverVectors;
    m_graverVectors = new GraVec (*m_projected_lattice_bases[m_rank], false);
}

void
GraverComputeState::MakeGraverBasisWithZSolve() 
{
    VectorArray v( get_vectors() );
    MakeGraverBasisWithZSolve (v);
    delete m_graverVectors;
    m_graverVectors = new GraVec ( std::move(v) );
};

void 
GraverComputeState::liftOneStep () {
    std::cout << "Lifting one step: " << m_graverVectors->get_size() << "->" << m_graverVectors->get_size()+1 << "\n";
    std::cout << "Current Vectors:\n";
    std::cout << m_graverVectors->get_vectors();
    std::cout << "Basis for lifting: \n";
    std::cout << *m_projected_lattice_bases[m_graverVectors->get_size()+1];
    m_graverVectors->lift (*m_projected_lattice_bases [m_graverVectors->get_size()+1]);
    std::cout << "Done lifting one step\n";
}

void print_start(){
    static int i = 0;
    std::cout << "Start mark no. " << i++ << "\n";
}

void print_ende(){
    static int i = 0;
    std::cout << "End mark no. " << i++ << "\n";
}


void 
GraverComputeState::liftGraverProperty () {
    IntegerType max_norm = m_graverVectors->get_max_norm ();
    std::cout << "Maximum norm on first components:" << max_norm << "\n";

    typedef _4ti2_zsolve_::NormPair< IntegerType > NP;
    typedef std::vector< NP > NPvector;
    NPvector jobs;
    if (m_graverVectors->has_vectors_with_norm(1))
	jobs.push_back( NP(1,1) );
    IntegerType current_norm = 0;
    while (current_norm < 2*max_norm) {
	current_norm++;
	std::cout << "Now doing norm: "<< current_norm << "\n";
	std::vector < std::future < VectorArray > > m_futures;
	// std::vector < VectorArray > m_futures;
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
	// std::cout << "Lowest norm jobs : ";
	// for (auto it = lowest_norm_jobs.begin(); it != lowest_norm_jobs.end(); it++) {
	//     std::cout << "(" << it->first << "," << it->second << "),";
	// }
	// std::cout << "\n";
	// std::cout << "Other jobs: ";
	// for (auto it = jobs.begin(); it != jobs.end(); it++) {
	//     std::cout << "(" << it->first << "," << it->second << "),";
	// }
	// std::cout << "\n";
	// if (lowest_norm_jobs.size() == 0) {
	//     std::cout << "No jobs of norm "<< current_norm << "\n";
	// }
	// // Do lowest norm jobs:
	for (auto it = lowest_norm_jobs.begin(); it != lowest_norm_jobs.end(); it++) {
	    std::cout << "Doing job :" << it->first << "," << it->second << "\n";
	    // Check if there are vectors with this norm:
	    if ( m_graverVectors->has_vectors_with_norm( it->first ) &&
		 m_graverVectors->has_vectors_with_norm( it->second )){
		// Synchronous version:
		// m_futures.push_back ( graverJob2 (m_graverVectors->get_filter(it->first), m_graverVectors->get_filter(it->second)));
		// Async version: 
		// Remarks: 
		// - the std::cref wrappers are needed to prevent copy on pass (see http://stackoverflow.com/questions/14851163/why-does-stdasync-copy-its-const-arguments)
		// - the launch directive can be modified to start more or less threads, having no directive lets the runtime code decicde
//		try {
		// std::cout << "Starting job asynchronously \n";
		std::future < VectorArray > fut = std::async(
		    std::launch::async, // This directive makes it launch a new thread for each job (not good if there are many!)
		    &GraverComputeState::graverJob2,
		    this,
		    std::cref (m_graverVectors->get_filter(it->first)),
		    std::cref (m_graverVectors->get_filter(it->second)));
		m_futures.push_back (std::move(fut));
//		}
//		catch (std::system_error e) {
//		    std::cout << "Can't create threads :(, exiting.\n";
//		    std::exit (1);
//		}
		// For debug purposes, wait for finish?
		// fut.wait();
	    } // if there are pairs of vectors in this job
	} // for over lowest_degree_jobs
	// Synchronize:
  	std::cout << "Waiting for sync ... ";
	std::cout.flush();
  	for (auto it = m_futures.begin(); it != m_futures.end(); ++it)
  	    it->wait();
  	std::cout << "Done.\n";

	// Storing new Graver elements: No newly created vector is
	// redundant, but there may be duplicates among the vectors
	// retrived from the different jobs.  To remove them, we build
	// a reduction tree only for new vectors.
	ReductionTree *R = new ReductionTree;
	for (auto it = m_futures.begin(); it != m_futures.end(); ++it) {
	    VectorArray res = it->get();
	    // VectorArray res = std::move(*it);
	    if (res.get_number() == 0) // nothing new
		continue;
	    std::cout << "Current norm: " << current_norm << "\n";
	    // check if maximum increased
	    if (max_norm < current_norm) {
		std::cout << "New norm bound : " << current_norm << "\n";
		max_norm = current_norm;
	    }
	    // std::cout << "I'm going to add new vectors. So far I got " << m_graverVectors->get_number() << std::endl;

	    std::cout << "Need to run " << res.get_number() << " reduction tests :(" << std::endl;
	    for (int j = 0; j<res.get_number(); j++){
		if (! R->isReducible(res[j])) {
		    R->insert (res[j]);
		    m_graverVectors->insert(std::move(res[j]));
		}
	    }
	    std::cout << "Done with that" << std::endl;
	    // std::cout << "and now there are: " << m_graverVectors->get_number() << std::endl;
	}
	// Clean up futures:
	m_futures.clear();
	delete R;
	    
	// Now all jobs with norm sum <= current_norm are done.
	// Add new jobs for each norm pair (i, current_norm), i=
	// 1..current_norm such that there are moves in the
	// respective degrees.
	// std::cout << "Current Jobs overview: \n";
// 	for (auto it=jobs.begin(); it != jobs.end(); it++) {
// 	    std::cout <<"("<<it->first<<","<<it->second<<"), ";
// 	}
	// std::cout << std::endl;
	// std::cout << "Current size of Graver basis: " << m_graverVectors->get_number() << "\n";
	if (m_graverVectors->has_vectors_with_norm(current_norm))
	    for (IntegerType i = 1; i <= current_norm; i++)
		if (m_graverVectors->has_vectors_with_norm(i))
		    jobs.push_back ( NP (i, current_norm));
	// Keep jobs sorted according to total norm
	// std::cout << "Sorting remaining jobs ... ";
	// std::cout.flush()
	std::sort(jobs.begin(), jobs.end());
	// std::cout << "done!\n";
    }; // while (current_norm < 2*max_norm)
}

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


VectorArray
GraverComputeState::graverJob (const VecVecP& Gr, const VecVecP& Gs) const
{
    std::cout << "I'm doing a job of size " << Gr.size() * Gs.size() << "...";
    // Step 1: Compute all sums f + g where f and g are sign
    // consistent with f \in Gr, g \in Gs.
    VectorArray result (0,Gr[0].v->get_size());
    for (uint i = 0; i < Gr.size(); i++ ) {
	// if Gr and Gs are the same, then we don't need all pairs:
	uint j = 0;
	if (&Gr == &Gs)
	    j = i+1;
	for (; j < Gs.size(); j++ ) {
	    // Check for sign consistency on the first n-1 components and anti-sign consistenty on last entry
	    if ( Gr[i].is_applicable (Gs[j])) 
	    {
		// Todo: 'Quickcheck'
		Vector sum = *Gr[i].v + *Gs[j].v;
		std::cout << "I decided to check the sum " << *Gr[i].v << " + " << *Gs[j].v << "=" << sum << "\n" ;
		if (m_graverVectors->is_reducible (sum)) {
		    std::cout << "oops, was reducible... \n";
		}
		else {
		    std::cout << "great, not reducible ! \n";
		    result.insert(std::move(sum));
		}
		// result->insert (new Vector (Gr[i] + Gs[j]));
	    }
	}
    }
    std::cout << " ... Done.\n";
    return result;
}

// void dumpPattern (BitSet bs){
//     std::cout << bs << "\n";
// }

VectorArray
GraverComputeState::graverJob2 (const GraverFilter& Gr, const GraverFilter& Gs) const
{
    VectorArray result (0,(Gr.begin()->second)[0].v->get_size());
    for (auto it = Gr.begin(); it != Gr.end(); it++ ) {
	auto matching = Gs.begin();
	if (&Gs == &Gr) {
	    // scroll forward to avoid duplicates
	    std::advance (matching, std::distance (Gr.begin(), it));
	}
	for (; matching != Gs.end(); matching++) {
	    if (!  ( BitSet::set_disjoint (matching->first.first , it->first.second ) &&
		     BitSet::set_disjoint (matching->first.second, it->first.first  ) 
		    ))
		continue;
	    for (uint i = 0; i < it->second.size(); i++) {
		for (uint j = 0; j < matching->second.size(); j++) {
		    // std::cout << "On: " << *it->second[i].v << " + " << *matching->second[j].v << "\n" ;
		    // Check for sign inconsistency on the last component 
		    if (it->second[i].last_entry() <= 0 && matching->second[j].last_entry() <= 0 )
			continue;
		    if (it->second[i].last_entry() >= 0 && matching->second[j].last_entry() >= 0 )
			continue;
		    // Todo: 'Quickcheck'
		    Vector sum = *it->second[i].v + *matching->second[j].v;
		    // std::cout << "I decided to check the sum " << *it->second[i].v << " + " << *matching->second[j].v << "=" << sum << "\n" ;
		    if (m_graverVectors->is_reducible (sum)) {
			// std::cout << "was reducible... \n";
		    }		    
		    else {
			// std::cout << "great, not reducible ! \n";
			result.insert(std::move(sum));
		    }
		}
	    }
	}
    }
    // std::cout << " ... Done.\n";
    return result;
}
