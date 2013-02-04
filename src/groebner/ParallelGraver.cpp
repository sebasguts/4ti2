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

#include <iostream>
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
//    int num_variables = feasible.get_dimension();

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

    print_vector (P);
    // Project and compute a Graver basis with zsolve.
    

    // Undo the permutation so that the users coordinates are
    // restored.
    basis.permute(P);
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
