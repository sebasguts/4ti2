#include <gtest/gtest.h>

#include "../ValueTree.hpp"
#include "../Lattice.hpp"
#include "../Vector.hpp"

using namespace _4ti2_zsolve_;

class ValueTreeTest : public ::testing::Test {
protected:
    virtual void SetUp() {
    // Initialize ValueTrees if needed
	L1 = new VectorArray<int> (2,8);
	L1->append_vector( create_vector<int> (8,1));
	L1->append_vector( create_vector<int> (8,2));
	T1 = new ValueTree<int> (L1);
    }
    
    // virtual void TearDown() 
    ValueTree<int> *T1;
    VectorArray <int>* L1;
};

TEST_F(ValueTreeTest, EmptyTree) {
    EXPECT_EQ(NULL, T1->zero);
    EXPECT_EQ(-1, T1->level);
}

TEST_F(ValueTreeTest, Insert) {
    int *v1 = create_vector<int> (8, 1);
    int *v2 = create_vector<int> (8, 2);
    T1->insert_vector(0, true);
    T1->insert_vector(1, true);
    T1->dump();
}
