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

#ifndef _4ti2_zsolve__UniqueVectorsHash_
#define _4ti2_zsolve__UniqueVectorsHash_

#include <iostream>
#include <unordered_set>

#include <zsolve/Vector.hpp>

// Implements a set of zsolve vectors using std::unordered_set

// see http://stackoverflow.com/questions/15728266/fast-ways-to-remove-duplicates-from-a-list-of-integer-vectors
// for some discussion on using hashing to extract unique integer vectors.

namespace _4ti2_zsolve_ 
{

// zsolve vectors don't save their length.  That's why we wrap them.
template <typename T>
struct WrappedVector {
    T *p;
    size_t length;
};

template <typename T>
struct VectorHash {
    size_t operator()(const WrappedVector<T>& v) const 
    {
	size_t result = 2166136261U;
	for (size_t i = 0; i<v.length; i++){
	    result = 127 * result + static_cast<size_t>( (v.p)[i] );
	}
	return result;
    }
};

// For gmp integers we just use the overflowed blocks for hashing
template <>
struct VectorHash <mpz_class> {
    size_t operator()(const WrappedVector<mpz_class>& v) const 
    {
	size_t result = 2166136261U;
	for (size_t i = 0; i<v.length; i++){
	    result = 127 * result + static_cast<size_t>( (v.p)[i].get_ui() );
	}
	return result;
    }
};

template <typename T>
struct VectorEqual {
    bool operator()(const WrappedVector<T>& v1, const WrappedVector<T>& v2) const {
	for (size_t i=0; i<v1.length; i++){
	    if ((v1.p)[i] != (v2.p)[i])
		return false;
	}
	return true;
    }
};

template <typename T> 
class UniqueVectorsHash {

public:    
    UniqueVectorsHash (size_t _l) : length(_l) { };
    ~UniqueVectorsHash () { 
	for (auto it = hash_table.begin(); it != hash_table.end(); it++) {
	    delete_vector<T> ((*it).p);
	}
    }

private:
    std::unordered_set<WrappedVector<T>, VectorHash<T>, VectorEqual<T> > hash_table;
    size_t length;

public:
    void insert (T *v) {
	WrappedVector<T> w;
	w.p = v;
	w.length = length;
	hash_table.insert (w);
    }

    void insert_copy (T *v) {
	WrappedVector<T> w;
	w.p = copy_vector<T> (v, length);
	w.length = length;
	hash_table.insert (w);
    }

    bool is_present (T *v) {
	WrappedVector<T> w;
	w.p = v;
	w.length = length;
	// std::cout << "Hash value: " << VectorHash<T>()(w);
	return (hash_table.count(w) > 0);
    }

    bool is_present (WrappedVector<T> *v) {
	return (hash_table.count(v) > 0);
    }

    void dump() {
	std::cout << "Dump: \n";
	for (auto it = hash_table.begin(); it != hash_table.end(); it++){
	    print_vector (std::cout, it->p, length);
	    std::cout << "\n";
	}
	std::cout << "\n";
    }
};


}

#endif
