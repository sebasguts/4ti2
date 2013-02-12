#include <gtest/gtest.h>
#include <iostream>

typedef int64_t IntegerType;

#include "groebner/Feasible.h"
#include "groebner/FeasibleStream.h"
#include "groebner/LatticeBasis.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/ParallelGraver.h"
#include "groebner/Permutation.h"

using namespace _4ti2_;

class ParallelGraverTest : public ::testing::Test {
protected:
    virtual void SetUp() {
	v1 = new VectorArray(3,4,0);
	(*v1)[0][0] = 1;
	(*v1)[0][1] = 1;
	(*v1)[1][2] = 1;
	(*v1)[2][3] = 1;

	// content of A.mat
	// 1 4
	// 1 2 3 4
	Feasible* feasible = input_Feasible("A");
	v2 = new VectorArray(0, 4);
	lattice_basis (feasible->get_matrix(), *v2);
    }

    // virtual void TearDown() {}
    VectorArray *v1;
    VectorArray *v2;
};

TEST_F(ParallelGraverTest, readFile) {
    EXPECT_EQ(v2->get_number(), 3);
    EXPECT_EQ(v2->get_size(), 4);
}

TEST_F(ParallelGraverTest, projectAndLift) {
    VectorArray *my_v2 = new VectorArray (*v2); // copy
    ParallelGraver::permute_full_rank_to_left (*my_v2);
    VectorArray *my_v2_proj = new VectorArray (3, 3);
    // The projected lattice will be a 3x3 matrix
    VectorArray::project(*my_v2, 0,2, *my_v2_proj); // project
    delete my_v2;
    my_v2 = ParallelGraver::liftToIndex(*my_v2_proj, *v2, 3);
    EXPECT_EQ(*my_v2, *v2);
    delete my_v2_proj;
    delete my_v2;
}

TEST_F(ParallelGraverTest, permuteRank) {
    VectorArray *my_va = new VectorArray (*v1); //copy
    Permutation p = ParallelGraver::permute_full_rank_to_left (*my_va);
    EXPECT_FALSE(*my_va == *v1);
    my_va->permute(p);
    EXPECT_TRUE(*my_va == *v1);
    delete my_va;
}

TEST_F(ParallelGraverTest, signConsistent) {
    EXPECT_TRUE(ParallelGraver::sign_consistent ((*v2)[0], (*v2)[1], 4));
    Vector *u1 = new Vector (3);
    Vector *u2 = new Vector (3);
    (*u1)[0] = -1;
    (*u1)[1] = 199;
    (*u1)[2] = 3;
    (*u2)[0] = -4;
    (*u2)[1] = 2;
    (*u2)[2] = -3;
    EXPECT_TRUE(ParallelGraver::sign_consistent ((*u1), (*u2), 2));
    EXPECT_FALSE(ParallelGraver::sign_consistent ((*u1), (*u2), 3));
    delete u1;
    delete u2;
}

TEST(ParallelGraver, isBelow) {
    Vector *u1 = new Vector (3);
    Vector *u2 = new Vector (3);
    Vector *u3 = new Vector (3);
    Vector *u4 = new Vector (3, 0);
    (*u1)[0] = -1;
    (*u1)[1] = 199;
    (*u1)[2] = 3;
    (*u2)[0] = -4;
    (*u2)[1] = 2;
    (*u2)[2] = -3;
    (*u3)[0] = -3;
    (*u3)[1] = 2;
    (*u3)[2] = -3;
    EXPECT_FALSE(ParallelGraver::is_below(*u1, *u2));
    EXPECT_FALSE(ParallelGraver::is_below(*u2, *u1));
    EXPECT_TRUE(ParallelGraver::is_below(*u3,*u2));
    EXPECT_TRUE(ParallelGraver::is_below(*u4,*u1));
    EXPECT_TRUE(ParallelGraver::is_below(*u4,*u2));
    EXPECT_TRUE(ParallelGraver::is_below(*u4,*u3));
    delete u1;
    delete u2;
}
