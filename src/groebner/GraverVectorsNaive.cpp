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

#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/GraverAPI.hpp"

#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/GraverVectors.h"
#include "groebner/GraverVectorsNaive.h"

using namespace _4ti2_;

GraverVectorsNaive::GraverVectorsNaive () {
}

GraverVectorsNaive::GraverVectorsNaive (const VectorArray& va){
    m_data = new VectorArray (va);
};
    
GraverVectorsNaive::~GraverVectorsNaive() {
    delete m_data;
}

GraverVectorsIterator
GraverVectorsNaive::begin (IntegerType norm){
    // Fix me
    return m_data2.begin();
};

GraverVectorsIterator
GraverVectorsNaive::end (IntegerType norm){
    // Fix me
    return m_data2.end();
}

void 
GraverVectorsNaive::addNegatives () {
    // Fix me
}  

void 
GraverVectorsNaive::removeNegatives (bool lexicographic) {
    // Fix me
}

Size 
GraverVectorsNaive::get_size(){
    return m_data->get_size();
} 
    
void 
GraverVectorsNaive::insert (Vector v) {
    m_data->insert(v);
} 

void 
GraverVectorsNaive::insert (Vector&& v) {
    m_data->insert(std::move(v));
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
GraverVectorsNaive::lift_with_basis( 
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

void
GraverVectorsNaive:: lift (const VectorArray& lifted_basis)
{
    // This is a trick.  Since zsolve solves Ax = b, we need to
    // transpose the basis to get a coefficient vector (remember,
    // members of the basis are rows originally.)  We don't need to
    // transpose the lifted basis because lift_with_basis correctly
    // combines its rows according to the coefficients.
    VectorArray *basis = new VectorArray (lifted_basis.get_number(), lifted_basis.get_size()-1);
    VectorArray::project(lifted_basis, 0, lifted_basis.get_size()-1, *basis);
    VectorArray *basis_transposed = new VectorArray (lifted_basis.get_size(), basis->get_number());
    VectorArray::transpose (*basis, *basis_transposed);
    delete basis;

    VectorArray *result = new VectorArray (0, lifted_basis.get_size());
    for (int i = 0; i < m_data->get_number(); i++)
	// The following does a move since the return value of a
	// function binds as an rvalue reference (hopefully).
	result->insert( lift_with_basis ((*m_data)[i], *basis_transposed, lifted_basis) );
    delete basis_transposed;
    delete m_data;
    m_data = result;
}



