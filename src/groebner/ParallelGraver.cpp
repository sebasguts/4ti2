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
#include "groebner/LatticeBasis.h"
#include "groebner/HermiteAlgorithm.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/GraverComputeState.h"

#include <algorithm>
#include <iostream>

using namespace _4ti2_;

ParallelGraver::ParallelGraver()
{
}

ParallelGraver::~ParallelGraver()
{
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
    // Dimension of the problem
    int num_variables = feasible.get_dimension();
    *out << "Number of variables of problem: " << num_variables << "\n";

    *out << "Computing Graver basis (parallel) ...\n";

    // Store the lattice basis
    lattice_basis (feasible.get_matrix(), basis);

    // Rank of the lattice
    int m_rank = basis.get_number();  // Human numbering
    *out << "Rank of problem: " << m_rank << "\n";

    // Permute columns such that a full rank matrix is in the first
    // rank many columns.  P stores the permutation which undoes the
    // shuffling in the end.
    Permutation P = permute_full_rank_to_left (basis);

    // Create GraverComputeState.
    GraverComputeState current (basis);  // give lattice basis for initialization

    // Compute a Graver basis of the smallest projection with zsolve.
    current.projectToRank();
    current.MakeGraverBasisWithZSolve (); // Should also add in the negatives!
    
    *out << "Graver basis of the projected lattice: \n";
    *out << current.get_vectors();
    *out << "\n";

    // The big lifting loop. varindex stores the index of the next
    // variable in machine counting.  For instance if we have a kx3
    // matrix of already lifted stuff, then varindex would start out
    // as "3" and all norms to be computed in the following are over
    // 0,1,2.
    for (Index varindex = m_rank; varindex < num_variables; varindex++ ) {
	*out << "Now lifting variable number " << varindex+1 << "\n";

	current.liftOneStep ();

 	*out << "Graver basis property on rank many components:\n";
 	*out << current.get_vectors() << std::endl;
 	*out << "\n";

	current.liftGraverProperty ();
    } // End of the big variable lifting loop

    // Undo the permutation
    basis = current.get_vectors_without_negatives();
    basis.permute(P);
    
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


#if 0
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
	IntegerType current_norm = 0;
	while (current_norm < 2*max_norm) {
	    current_norm++;
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
		    // Remarks: 
		    // - the std::cref wrappers are needed to prevent copy on pass (see http://stackoverflow.com/questions/14851163/why-does-stdasync-copy-its-const-arguments)
		    // - the launch directive can be modified to start more or less threads, having no directive lets the runtime code decicde
		    try {
			std::future < VectorArray* > fut = std::async(
			    // std::launch::async, // This directive makes it launch a new thread for each job (not good if there are many!)
			    GraverComputeState::graverJob,
			    std::cref (*(current.get_tree().at(it->first))),
			    std::cref (*(current.get_tree().at(it->second))),
			    std::cref (current.get_vectors()),
			    varindex
			    );
			m_futures.push_back (std::move(fut));
		    }
		    catch (std::system_error e) {
			std::cout << "Can't create threads :(, exiting.\n";
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
	    
	    // Now all jobs with norm sum <= current_norm are done.
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
	// std::cout << "Ok buddy, I'm done with this variable, here's what I got\n";
	// std::cout << current.get_vectors();
	// TODO: At this point "current_graver_basis" will be destructed because the destructor of the 
#endif
