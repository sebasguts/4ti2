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
#include "groebner/VectorArrayStream.h"

namespace _4ti2_
{

// GraverComputeState::GraverComputeState (const VectorArray& vs){
//     m_normTree = new NormBST (); // Does not actually compute the tree
//     m_normArray = new std::vector<IntegerType>;	
//     m_array = new VectorArray (vs);
// }

// TODO: Implement proper move semantics (need support in VectorArray)
// GraverComputeState::GraverComputeState (VectorArray&& vs){
//     m_normTree = new NormBST ();
//     m_normArray = new std::vector<IntegerType>;	
//     m_array = new VectorArray (std::move(vs)); 
// }

GraverComputeState::GraverComputeState (VectorArray *vs){
    m_normTree = new NormBST ();
    m_latticeBasis = new VectorArray (0,vs->get_size());
    m_array = vs; // No ownership !
}

GraverComputeState::~GraverComputeState() {
    // do _NOT_ delete m_array, it is not yours!
    delete m_normTree;
    delete m_latticeBasis;
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
	// The map container gives a strong guarentee of not being
	// changed if an exception occurs
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
}

IntegerType 
GraverComputeState::maximum_norm () const {
    assert(m_normTree != NULL);
    return m_normTree->rbegin()->first;
}

} // namespace _4ti2_
