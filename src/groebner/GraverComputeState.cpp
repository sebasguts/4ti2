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
#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"

#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/GraverAPI.hpp"

namespace _4ti2_
{

GraverComputeState::GraverComputeState (const VectorArray& lb){
    m_normTree = new NormBST (); // Does not actually compute the tree
    m_array = new VectorArray (0, lb.get_size());
    m_latticeBasis = new VectorArray (lb);

    for (int i = lb.get_number(); i < lb.get_size(); i++){
	m_projected_lattice_bases.push_back (new VectorArray(lb.get_number(), i));
	VectorArray::project(lb, 0, i, *m_projected_lattice_bases.back());
	// std::cout << "Projected basis: " << *m_projected_lattice_bases.back();
    }
    // Insert the full basis too for lifting convenience
    m_projected_lattice_bases.push_back(new VectorArray (lb));
}

// GraverComputeState::GraverComputeState (VectorArray *vs){
//     m_normTree = new NormBST ();
//     m_latticeBasis = new VectorArray (0,vs->get_size());
//     m_array = vs; // No ownership !
// }

GraverComputeState::~GraverComputeState() {
    delete m_array;
    delete m_normTree;
    delete m_latticeBasis;
    for (auto it = m_projected_lattice_bases.begin(); it != m_projected_lattice_bases.end(); it++)
	delete *it;
}

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

/** 
 * \brief Lift a Vector according to given lift on bases.
 * 
 * @param v Vector to be lifted
 * @param basis VectorArray consisting of the projections of lifted_basis
 * @param lifted_basis VectorArray of long basis vectors, aligned with basis
 *
 * see http://stackoverflow.com/questions/4986673/c11-rvalues-and-move-semantics-confusion
 * for why we (seem to) return by value.
 * 
 * @return A new Vector that is the lift
 */
Vector
GraverComputeState::lift_with_basis( 
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
    
    Vector result_vector (lifted_basis.get_size(), 0);
    for (int i = 0; i < lifted_basis.get_number(); i++)
	result_vector.add( (lifted_basis[i]) * ((*coefficient_vector)[i]));
    delete coefficient_vector;
    return result_vector;
}

VectorArray
GraverComputeState::lift_with_basis( 
    const VectorArray& va,
    const VectorArray& basis,
    const VectorArray& lifted_basis)
{
    // This is a trick.  Since zsolve solves Ax = b, we need to
    // transpose the basis to get a coefficient vector (remember,
    // members of the basis are rows originally.)  We don't need to
    // transpose the lifted basis because lift_with_basis correctly
    // combines its rows according to the coefficients.
    VectorArray *basis_transposed = new VectorArray (basis.get_size(), basis.get_number());
    VectorArray::transpose (basis, *basis_transposed);

    VectorArray result = VectorArray (0, lifted_basis.get_size());
    for (int i = 0; i < va.get_number(); i++)
	// The following does a move since the return value of a
	// function binds as an rvalue reference (hopefully).
	result.insert(lift_with_basis (va[i], *basis_transposed, lifted_basis));
    delete basis_transposed;
    return result;
}

void
GraverComputeState::liftToIndex ( Index index)
{
    // Interesting... does this leak memory?  It should not since
    // VectorArray's move assignment swaps the pointers and thus the
    // old m_array gets cleaned up.
    *m_array = lift_with_basis (*m_array,
				*m_projected_lattice_bases[index],
				*m_projected_lattice_bases[index+1]);
}


IntegerType 
GraverComputeState::maximum_norm () const {
    assert(m_normTree != NULL);
    // std::cout << m_normTree->rbegin()->first << "\n";
    return m_normTree->rbegin()->first;
}

} // namespace _4ti2_
