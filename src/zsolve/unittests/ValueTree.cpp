#include <gtest/gtest.h>

#include "../ValueTree.hpp"

class ValueTreeTest : public ::testing::Test {
protected:
    virtual void SetUp() {
	// Initialize ValueTrees if needed
    }

    // virtual void TearDown() {}

    _4ti2_zsolve_::ValueTree<int> t0_;
    _4ti2_zsolve_::ValueTree<int> t1_;
    _4ti2_zsolve_::ValueTree<int> t2_;
};

TEST_F(ValueTreeTest, EmptyTree) {
    EXPECT_EQ(NULL, t0_.zero);
    EXPECT_EQ(-1, t0_.level);
}
