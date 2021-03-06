/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Peter Malkin.

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

#ifndef _4ti2_qsolve__BinaryTree_
#define _4ti2_qsolve__BinaryTree_

#include <vector>

namespace _4ti2_
{

template <class IndexSet>
class BinaryTree
{
public:
    BinaryTree();
    ~BinaryTree();

    void insert(const IndexSet& support, Index index);
    void insert(const std::vector<IndexSet>& supports);
    void insert(const std::vector<IndexSet>& supports, Index start, Index end);
    bool dominated(const IndexSet& b, Index index1, Index index2) const;
    void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
    void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
    void find(std::vector<std::pair<Index,Index> >& inds, const BinaryTree& tree, Size count) const;
    void clear();
    void dump() const;

private:
    class TreeNode
    {
    public:
        TreeNode() {}
        virtual ~TreeNode() {}
        virtual TreeNode* insert(const IndexSet& support, Index index) = 0;
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const = 0;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const = 0;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const = 0;
        virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, const TreeNode* tree, Size count) const = 0;
        virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, Index index, Size count) const = 0;
        virtual void dump(int level) const = 0;
    };

#if 0
    class TreeLeaf : public TreeNode {
    public:
        TreeLeaf(Index _index, const IndexSet& _is) : i(_index), is(_is) {}
        ~TreeLeaf() {}
        TreeNode* insert(const IndexSet& support, Index index);
        bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const { inds.push_back(i); }
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const { inds.push_back(i); }
        //virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, const TreeNode* tree, Size count) const;
        //virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, Index index, Size count) const;
        void dump(int level) const;
    private:
        Index i;
        IndexSet is;
    };
#endif

    class TreeBranch : public TreeNode {
    public:
        TreeBranch(Index index, const IndexSet& is0, TreeNode* _zero, const IndexSet& is1, TreeNode* _one);
        ~TreeBranch();
        TreeNode* insert(const IndexSet& support, Index index);
        bool dominated(const IndexSet& s, Index index1, Index index2) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, const TreeNode* tree, Size count) const;
        virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, Index index, Size count) const;
        void dump(int level) const;
    private:
        Index i;
        IndexSet is0;
        IndexSet is1;
        TreeNode* zero;
        TreeNode* one;
    };

#if 1
    class TreeBucket: public TreeNode {
    public:
        TreeBucket();
        ~TreeBucket();
        TreeNode* insert(const IndexSet& support, Index index);
        bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, const TreeNode* tree, Size count) const;
        virtual void find(std::vector<std::pair<Index,Index> >& inds, const IndexSet& supp, Index index, Size count) const;
        void dump(int level) const;
    private:
        std::vector<IndexSet> supps;
        std::vector<Index> indices;
        static const Size MAX_BUCKET_SIZE;
    };
#endif
    static const int INDENT;

    TreeNode* root;
    Size supp_size;
};

} // namespace _4ti2_

// Include template defitinions.
#include "qsolve/BinaryTree.hpp"

#endif
