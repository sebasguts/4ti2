#include <gtest/gtest.h>
#include <iostream>

typedef int64_t IntegerType;

#include "groebner/ReductionTree.h"

using namespace _4ti2_;

class ReductionTreeTest : public ::testing::Test {
protected:
    virtual void SetUp() {
	v = Vector (5, 1);
	w = Vector (5, 2);
	u = Vector (5, -3);
	x = Vector (5, 0);
    }

    // virtual void TearDown() {}
    ReductionTree T1, T2;
    Vector u, v, w, x;
};


TEST_F(ReductionTreeTest, insert) {
    T1.insert(u);
    T1.insert(v);
    T1.insert(w);
    T1.insert(x);
}

TEST_F(ReductionTreeTest, basic_reduction) {
    T1.insert(v);
    EXPECT_FALSE(T1.isReducible (u));
    EXPECT_TRUE(T1.isReducible (v));
    v[0] = 0;
    EXPECT_FALSE(T1.isReducible(v));
    T1.insert(u);
    u[0] = -2;
    EXPECT_FALSE(T1.isReducible(u));
    u[0] = -4;
    EXPECT_TRUE(T1.isReducible(u));
    u[1] = 1;
    EXPECT_FALSE(T1.isReducible(u));
}

TEST_F(ReductionTreeTest, clear) {
    T1.insert(v);
    EXPECT_TRUE(T1.isReducible(v));
    T1.clear();
    EXPECT_FALSE(T1.isReducible(v));
}

TEST_F(ReductionTreeTest, zero_in_reduction) {
    T1.insert(x); // insert the zero-vector which means that everything is reducible
    EXPECT_TRUE(T1.isReducible(u));
    EXPECT_TRUE(T1.isReducible(v));
    EXPECT_TRUE(T1.isReducible(w));
    T2.insert (u);
    T2.insert (v);
    T2.insert (w);
    // The zero vector can not be reduced by non-zero vectors
    EXPECT_FALSE(T2.isReducible(x));
}

TEST_F(ReductionTreeTest, fromBugs) {
    Vector v1 (3,0);
    Vector v2 (3,0);
    Vector v3 (3,0);
    Vector v4 (3,0);
    v1[0]=1; v1[1]= 1; v1[2]= -1 ;
    v1[0]=1; v2[1]=-1; v2[2]= -1 ;
    v1[0]=2; v3[1]= 0; v3[2]= -2 ;
    v1[0]=0; v4[1]= 2; v4[2]=  0 ;
    T1.insert(v1);
    T1.insert(v2);
    T1.insert(v3);
    T1.insert(v4);
    Vector t (3,0);
    t[0]=0;t[1]=3;t[2]=0;
    EXPECT_TRUE(T1.isReducible(t));
}
