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

#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/GraverVectors.h"
#include "groebner/GraverVectorsWithFilter.h"
#include "groebner/BitSet.h"
#include "groebner/BitSetStream.h"
#include "groebner/LatticeBasis.h"

using namespace _4ti2_;

void
insertIntoFilter (GraverFilter& filter, const GraverVector& g){
    auto loc = filter.find ( std::pair<BitSet, BitSet> (g.pos, g.neg));
    if ( loc == filter.end() ){
	VecVecP tmp;
	tmp.push_back (g); // GraverVectors are cheap to copy
	std::pair<BitSet, BitSet> supports (g.pos, g.neg);
	filter.insert (
	    std::pair <  std::pair <BitSet, BitSet>,  VecVecP > (supports, std::move(tmp)));
    }
    else {
	loc->second.push_back(g);
    }
}

GraverVectorsWithFilter::GraverVectorsWithFilter () {
}

GraverVectorsWithFilter::GraverVectorsWithFilter (const VectorArray& va){
    m_data = new VectorArray (va);
    // size = m_data->get_size();
    createNormFilter();
    createRedTree();
};
    
GraverVectorsWithFilter::~GraverVectorsWithFilter() {
    delete m_data;
}

void 
GraverVectorsWithFilter::addNegatives () {
    // Fix me
}  

void 
GraverVectorsWithFilter::removeNegatives (bool lexicographic) {
    // Fix me
}

Size 
GraverVectorsWithFilter::get_size(){
    return m_data->get_size();
}

void printMark () {
    static int i = 1;
    std::cout << "Mark no: " << i++ << "\n";
}

void 
GraverVectorsWithFilter::insert (Vector&& v) {
    IntegerType norm = v.norm(m_data->get_size()-1);
    std::cout << "Running insert on norm " << norm << "\n";
    m_redTree.insert (v);
    m_data->insert(std::move (v));
    GraverVector g ( & (*m_data)[m_data->get_number()-1] );
    auto it = m_normFilter.find(norm);
    if (it == m_normFilter.end()) {
	// If the norm does not exist, create it
	// std::pair <IntegerType, std::vector< Vector* > > f(current_norm, );
	GraverFilter tmp;
	insertIntoFilter (tmp, g);
	m_normFilter.insert(std::pair <IntegerType, GraverFilter> (norm, std::move(tmp)));
    }
    else {
	// found
	insertIntoFilter( it->second, g);
    }
}

void 
GraverVectorsWithFilter::insert (VectorArray&& va) {
    std::cout << "Inserting an array ... \n";
    for (int i = 0; i < va.get_number(); i++){
	insert (std::move(va[i]));
    }
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
GraverVectorsWithFilter::lift_with_basis( 
    const Vector& v,
    const VectorArray& basis,
    const VectorArray& lifted_basis)
{
    // need to solve (basis * x = v)
    Vector sol (basis.get_number());
    solve (basis, v, sol);
    Vector result_vector (lifted_basis.get_size(), 0);
    for (int i = 0; i < lifted_basis.get_number(); i++)
	result_vector.add( (lifted_basis[i]) * (sol[i]));
    return result_vector;
}

void
GraverVectorsWithFilter::lift (const VectorArray& lifted_basis)
{
    // This is a trick.  Since zsolve solves Ax = b, we need to
    // transpose the basis to get a coefficient vector (remember,
    // members of the basis are rows originally.)  We don't need to
    // transpose the lifted basis because lift_with_basis correctly
    // combines its rows according to the coefficients.
    VectorArray *basis = new VectorArray (lifted_basis.get_number(), lifted_basis.get_size()-1);
    VectorArray::project(lifted_basis, 0, lifted_basis.get_size()-1, *basis);
    VectorArray *basis_transposed = new VectorArray (lifted_basis.get_size()-1, basis->get_number());
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
    // Todo: Do better here, those structures can be lifted too!
    createNormFilter();
    createRedTree(); 
}

void
GraverVectorsWithFilter::createRedTree() {
    m_redTree.clear();
    for (int i = 0; i < m_data->get_number(); i++)
	m_redTree.insert( (*m_data)[i] );
}

void dumpGraverFilter (const GraverFilter& gf) {
    for (auto it = gf.begin(); it != gf.end(); ++it){
	auto p = it->first;
	std::cout << "Index sets: \n" << p.first << " -- " << p.second << ": ";
	std::cout << it->second.size() << " vectors\n";
    }
}

void 
GraverVectorsWithFilter::dumpFilters (){
    for (auto it = m_normFilter.begin(); it != m_normFilter.end(); ++it){
	std::cout << "Dumping filter for norm " << it->first << "\n";
	dumpGraverFilter(it->second);
	std::cout << "---------------------\n";
    }
}

void
GraverVectorsWithFilter::createNormFilter () {
    std::cout << "Creating new norm filter" << std::endl;
    m_normFilter.clear();
    /// @TODO: Parallelize this computation by splitting the todolist
    for (int i = 0; i < m_data->get_number(); i++){
	std::cout << (*m_data)[i] << ": ";
	IntegerType current_norm = (*m_data)[i].norm(m_data->get_size()-1); /// Compute norm of first n-1 entries!
	std::cout << current_norm << "\n";
	GraverVector g ( &(*m_data)[i] );
	auto it = m_normFilter.find(current_norm);
	if (it == m_normFilter.end()) {
	    // If the norm does not exist, create it
	    // std::pair <IntegerType, std::vector< Vector* > > f(current_norm, );
	    GraverFilter tmp;
	    insertIntoFilter (tmp, g);
	    m_normFilter.insert(std::pair <IntegerType, GraverFilter> (current_norm, std::move(tmp)));
	}
	else {
	    // found
	    insertIntoFilter( it->second, g);
	}
    }
    std::cout << "Norm filter created, minimum norm: " << m_normFilter.begin()->first;
    std::cout << ", maximum norm : " << m_normFilter.rbegin()->first << "\n";
    dumpFilters();
}
