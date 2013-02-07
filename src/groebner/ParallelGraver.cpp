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

    // Project and compute a Graver basis with zsolve.
    VectorArray *projected = new VectorArray (basis.get_number(), m_rank);
    VectorArray::project(basis, 0, m_rank-1, *projected);

    *out << "Projected Lattice\n";
    *out << *projected << std::endl;

    // Call old zsolve Graver code to determine minimal elements.  We
    // need a lattice basis in zsolve Format: 
    
////     // @TODO Make precision
////     //customizable
////     _4ti2_zsolve_::GraverAPI <IntegerType> 
//// 	*m_graver_API = new _4ti2_zsolve_::GraverAPI <IntegerType>;
////     // The following is really stupid, but I believe it is a design
////     // error in the API.  I feel like I'm intended to only use
////     // _4ti2_matrix* as the return type of create_matrix.  However,
////     // that interface has no means of doing anything useful with the
////     // matrix.  In fact, if we look at the code of create_matrix then
////     // we see that a VectorArrayAPI is created and then upcasted to
////     // _4ti2_matrix.  For now we just cast it back.
////     _4ti2_zsolve_::VectorArrayAPI <IntegerType> *projected_g_basis_API = 
//// 	(_4ti2_zsolve_::VectorArrayAPI <IntegerType>*) m_graver_API->create_matrix(projected->get_number(), projected->get_size(), "lat");
////     _4ti2_zsolve_::VectorArray <IntegerType>& projected_g_basis = projected_g_basis_API->data;
////     for (int i = 0; i < projected->get_number(); i++) {
//// 	// Create a zsolve copy of this vector:
//// 	IntegerType *v = _4ti2_zsolve_::create_vector<IntegerType> (projected->get_size());
//// 	for (int j = 0; j < projected->get_size(); j++)
//// 	    v[j] = (*projected)[i][j];
//// 	projected_g_basis.append_vector(v);
////     }
////     // Now we have the projected lattice basis in a zsolve
////     // VectorArray, stored in the Graver_API.
////     m_graver_API->compute();
//// 
////     // Result is stored in zhom_data
////     
////     // Clean-up
////     delete projected;
////     delete m_graver_API;

    _4ti2_state *m_state = new _4ti2_zsolve_::GraverAPI <IntegerType> ();
    // _4ti2_state_set_options(m_state, argv, argc);
    _4ti2_matrix *m_matrix = m_state->create_matrix(projected->get_number(),
						    projected->get_size(),
						    "lat");
    // Fill matrix using the GMP API
    /// @TODO This will break if gmp is not available?
    for (int i = 0; i < projected->get_number(); i++) {
  	for (int j = 0; j < projected->get_size(); j++)
  	    m_matrix->set_entry_mpz_class (i,j, (*projected)[i][j]);
    }
    *out << "Matrix: \n";
    for (int i = 0; i < m_matrix->get_num_rows(); i++) {
	for (int j = 0; j < m_matrix->get_num_cols(); j++) {
	    int64_t temp_val = 0; 
	    m_matrix->get_entry_int64_t(i,j, temp_val);
	    *out << temp_val;
	}
	*out << std::endl;
    }
    m_state->compute();
    _4ti2_matrix *m_result = m_state->get_matrix("zhom");

    for (int i = 0; i < m_result->get_num_rows(); i++) {
	for (int j = 0; j < m_result->get_num_cols(); j++) {
	    int64_t temp_val = 0; 
	    m_result->get_entry_int64_t(i,j, temp_val);
	    *out << temp_val;
	}
	*out << std::endl;
    }
    delete m_state;

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
