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

#include "groebner/ReductionTree.h"

using namespace _4ti2_;

ReductionTree::ReductionTree()
{
}

ReductionTree::~ReductionTree() {
}

void
ReductionTree::clear () {
    m_branches.clear();
};

void 
ReductionTree::insert (const Vector& v) {
    insert_with_offset(v, 0);
}

void
ReductionTree::insert_with_offset (const Vector& v, const int offset) {
    if (offset == v.get_size()-1) {
	m_branches.insert ( std::pair <IntegerType, ReductionTree> (v[offset], ReductionTree() ));
    }
    else {
	auto loc = m_branches.insert (std::pair <IntegerType, ReductionTree> (v[offset], ReductionTree() ));
	// From the manual: insert returns a pair, with its member
	// pair::first set to an iterator pointing to either the newly
	// inserted element or to the element with an equivalent key
	// in the map. The pair::second element in the pair is set to
	// true if a new element was inserted or false if an
	// equivalent key already existed.
	loc.first->second.insert_with_offset(v, offset+1);
    }
}

bool 
ReductionTree::isReducible (const Vector& v) const {
    return isReducible_with_offset (v, 0);
}

bool 
ReductionTree::isReducible_with_offset (const Vector& v, const int offset) const {
    // std::cout << "in reduction with offset " << offset << " and component " << v[0] << std::endl;
    if (v[offset] == 0) {
	// we can only recurse on the zero branch:
	auto loc = m_branches.find( 0 );
	if (loc == m_branches.end())
	    return false;
	else 
	    return loc->second.isReducible_with_offset (v, offset+1);
    }
     else { // v[offset] != 0
 	IntegerType lower = v[offset] < 0 ? v[offset] : 0;
 	IntegerType upper = v[offset] > 0 ? v[offset] : 0;
 	// We find the first element that is larger than or equal to the
 	// lower bound, and then ask if it is smaller than or equal to the
 	// upper bound:
 	// std::cout << "Available branches: ";
 	// for (auto it = m_branches.begin(); it != m_branches.end(); it++) {
 	//     std::cout << it->first << " ";
 	// }
 	// std::cout <<"\n";
 	auto it = m_branches.lower_bound(lower);
 	// If it points to the end now, then we are done anyway:
 	if (it == m_branches.end()) {
 	    // std::cout << "Nothing beyond " << lower << "found. Not reducible\n";
 	    return false;
 	}
 	if (offset == v.get_size()-1) {
 	    if (it->first > upper) 
 		return false;
 	    else 
 		return true;
 	}
 	else {
 	    // recurse on each tree in the range between upper and lower
 	    for (; it != m_branches.upper_bound(upper); it++){
 		if (it->second.isReducible_with_offset(v, offset+1))
 		    return true;
 	    }
 	    // std::cout << "Key not found, breaking reduction with false.\n";
 	    return false;
 	}
     } // v[offset] != 0
}


