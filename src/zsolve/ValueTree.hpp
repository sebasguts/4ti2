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

namespace _4ti2_zsolve_
{

template <typename U> class ValueTree;

template <typename U> class ValueTreeNode
{
public:
    ValueTree <U> * sub_tree;
    U value;

    ValueTreeNode (U v, size_t vid) {
	sub_tree = new ValueTree <U> ();
	sub_tree->vector_indices.push_back (vid);
	value = v;
    }
    
    ~ValueTreeNode () {
	delete sub_tree;
    }
};


/**
 * \brief A binary tree with an integer value at each node
 *
 * A ValueTree stores integer vectors according to the norm of
 * projections.
 */
template <typename U> class ValueTree
{
public:
    int level;
    ValueTree <U>* zero;
    std::vector<ValueTreeNode <U> *> pos, neg;
    std::vector<size_t> vector_indices;
    
    ValueTree () {
	level = -1;
	zero = NULL;
    }

    ~ValueTree ()  {
	if (zero != NULL)
	    delete zero;
	for (size_t i = 0; i < pos.size (); i++)
	    delete pos[i];
	for (size_t i = 0; i < neg.size (); i++)
	    delete neg[i];
    }
};


} // namespace _4ti2_zsolve_

#endif
