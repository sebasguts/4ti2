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

    // Store the lattice generators
    basis = feasible.get_basis();
    // compute normal form
    int m_rank = hermite (basis, basis.get_size());
    *out << "Rank of problem: " << m_rank << "\n";
    // drop zero rows
    VectorArray temp (m_rank, basis.get_size());
    for (int i = 0; i < basis.get_number(); i++){
	if (! basis[i].is_zero())
	    temp[i] = std::move(basis[i]);
    }
    basis = std::move(temp);

    // The following would not work since the lattice may not be
    // saturated.  In this case get_matrix returns the matrix for the
    // saturation.  
    // lattice_basis (feasible.get_matrix(), basis);

    // Permute columns such that a full rank matrix is in the first
    // rank many columns.  P stores the permutation which undoes the
    // shuffling in the end.
    Permutation P = permute_full_rank_to_left (basis);
    // std::cout << "Permuted basis \n" << basis;

    // Create GraverComputeState.
    GraverComputeState current (basis);  // give lattice basis for initialization

    // Compute a Graver basis of the smallest projection with zsolve.
    current.projectToRank();
    // *out << current.get_vectors() << std::endl;
    current.MakeGraverBasisWithZSolve (); // Also adds the negatives!
    
//    *out << "Graver basis of the projected lattice: \n";
//    *out << current.get_vectors();
//    *out << "\n";

    // The big lifting loop. varindex stores the index of the next
    // variable in machine counting.  For instance if we have a kx3
    // matrix of already lifted stuff, then varindex would start out
    // as "3" and all norms to be computed in the following are over
    // 0,1,2.
    for (Index varindex = m_rank; varindex < num_variables; varindex++ ) {
	*out << "Now lifting variable number " << varindex+1 << "\n";

	current.liftOneStep ();

 	// *out << "Graver basis property on rank many components:\n";
  	// *out << current.get_vectors() << std::endl;
 	// *out << "\n";

	current.liftGraverProperty ();
    } // End of the big variable lifting loop

    // Undo the permutation
    std::cout << "removing negatives.... ";
    std::cout.flush();
    basis = current.get_vectors_without_negatives();
    std::cout << "done\n";
    basis.permute(P);
    
    std::cout << "all done\n";
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
#endif
