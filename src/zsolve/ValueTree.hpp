/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2013 4ti2 team.
Main author(s): Matthias Walter, Thomas Kahle.

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

#ifndef _4ti2_zsolve__ValueTree_
#define _4ti2_zsolve__ValueTree_

#include "zsolve/VectorArray.hpp"
#include "zsolve/Vector.hpp"

namespace _4ti2_zsolve_
{

/**
 * A value stores T vectors in a hierarchical way.  At level i it
 * stores all possible values that occur among the i-components of the
 * set of vectors.  It also stores a sub-tree for further down
 * components starting that is, the possible values for i+1
 * 
 */


template <typename T> class ValueTree;

template <typename T> class ValueTreeNode
{
public:
    ValueTree <T> * sub_tree;
    T value;

    ValueTreeNode (T v, size_t vid) {
	sub_tree = new ValueTree <T> ();
	sub_tree->vector_indices.push_back (vid);
	value = v;
    }
    
    ~ValueTreeNode () {
	delete sub_tree;
    }
};


/**
 * \brief A tree to store VectorArrays differently
 *
 * A ValueTree stores integer vectors according to the norm of
 * projections.
 */
template <typename T> class ValueTree
{

private:
    VectorArray <T> * myVectors; /// Pointer to set of vectors that this ValueTree stores.  Not owned by the valueTree!

public:
    int level; /// Which component does this value tree represent
    ValueTree <T>* zero;  /// Subtree for vectors that have an entry zero here
    std::vector<ValueTreeNode <T> *> pos, neg; /// Subtrees for positive and negative entries here
    std::vector<size_t> vector_indices; /// Auxilliary information: Which index did a vector have in some outside array.
    
    ValueTree () {
	level = -1;
	zero = NULL;
	myVectors = NULL;
    }

    ValueTree (VectorArray <T> *vectors) {
	level = -1;
	zero = NULL;
	myVectors = vectors;
    }

    ~ValueTree ()  {
	if (zero != NULL)
	    delete zero;
	for (size_t i = 0; i < pos.size (); i++)
	    delete pos[i];
	for (size_t i = 0; i < neg.size (); i++)
	    delete neg[i];
    }

public:
    void insert_vector (T* vector, size_t length, size_t vector_index, bool split_recursive);
    void insert_vector (size_t vector_index, bool split_recursive);
    void split (int current_variable, int start = -1);
    void dump ();
};

/// If the lattice is stored then we can infer some parameters:
template <typename T>
void 
ValueTree<T>::insert_vector (size_t vector_index, bool split_recursive) {
    assert (myVectors != NULL);
    insert_vector ((*myVectors)[vector_index], myVectors->num_variables(), vector_index, split_recursive);
}

template <typename T>
void 
ValueTree<T>::insert_vector (T* vector, size_t length, size_t vector_index, bool split_recursive) {
    if (level < 0) { 
	// Empty tree case 
	vector_indices.push_back (vector_index);
	if (split_recursive)
	    split();
    }
    else {
	T value = vector[level];
	if (value > 0)
	{
	    auto iter = pos.begin();
	    for (; iter != pos.end (); iter++)
		if (value <= (*iter)->value)
		    break;
	    if (iter == pos.end() || value != (*iter)->value)
		pos.insert (iter, new ValueTreeNode <T> (value, vector_index));
	    else
		(*iter)->sub_tree->insert_vector (vector, length, vector_index, split_recursive);
	}
	else if (value < 0)
	{
	    auto iter = neg.begin();
	    for (; iter != neg.end (); iter++)
		if (value >= (*iter)->value)
		    break;
	    if (iter == neg.end() || value != (*iter)->value)
		neg.insert (iter, new ValueTreeNode <T> (value, vector_index));
	    else
		(*iter)->sub_tree->insert_vector (vector, length, vector_index, split_recursive);
	}
	else // value == 0
	{
	    if (zero == NULL)
		zero = new ValueTree <T> ();
	    zero->insert_vector (vector, length, vector_index, split_recursive);
	}
    }
}

// /**
//  * What does this actually do?
//  * 
//  */
// template <typename T>
// void 
// ValueTree<T>::split (int current_variable, int start) {
//     int compo = start < 0 ? current_variable : start;
//     bool has_pos, has_neg, has_zero;
// 
//     if (level >= 0)
// 	return;
// 
//     for (; start < (int) current_variable; start++)
//     {
// 	compo = start < 0 ? current_variable : start;
// 	has_pos = has_neg = has_zero = false;
// 	for (size_t i = 0; i < tree->vector_indices.size (); i++) // vector_indices.size(): How many vectors are stored in this tree
// 	{
// 	    T value = (*m_lattice)[tree->vector_indices[i]][compo]; // Queries the current list of vectors :(
// 	    if (value > 0)
// 		has_pos = true;
// 	    else if (value < 0)
// 		has_neg = true;
// 	    else
// 		has_zero = true;
// 	    if (has_pos && has_neg)
// 		break;
// 	}
// 	if (has_pos && has_neg)
// 	    break;
//     }
//     if ((start < (int) current_variable) && (tree->vector_indices.size () >= 1))
//     {
// //            std::cout << "Splitting on " << compo << std::endl;
// 	tree->level = compo;
// 	for (size_t i = 0; i < tree->vector_indices.size (); i++)
// 	    insert_tree (tree, tree->vector_indices[i], false);
// 	start++;
// 	if (tree->zero != NULL)
// 	    split_tree (tree->zero, start);
// 	for (size_t i = 0; i < tree->pos.size (); i++)
// 	    split_tree (tree->pos[i]->sub_tree, start);
// 	for (size_t i = 0; i < tree->neg.size (); i++)
// 	    split_tree (tree->neg[i]->sub_tree, start);
//     }
// }


template <typename T>
void 
ValueTree<T>::dump ()
{
    if (level < 0)
    {
	std::cout << "dump::leaf:" << std::endl;
	for (size_t i = 0; i < vector_indices.size(); i++)
	{
	    std::cout << "  [" << vector_indices[i] << "] = ";
	    print_vector (std::cout, (*myVectors)[vector_indices[i]], myVectors->num_variables());
	    std::cout << std::endl;
	}
    }
    else
    {
	std::cout << "dump::node at level " << level << "\n";
	if (zero != NULL)
	{
	    std::cout << "dump::zero\n";
	    zero->dump();
	}
	for (size_t i = 0; i < pos.size (); i++)
	{
	    std::cout << "dump::pos (" << pos[i]->value << ")" << std::endl;
	    pos[i]->sub_tree->dump();
	}
	for (size_t i = 0; i < neg.size (); i++)
	{
	    std::cout << "dump::neg (" << neg[i]->value << ")" << std::endl;
	    neg[i]->sub_tree->dump();
	}
    }
}



} // namespace _4ti2_zsolve_

#endif
