#include <gtest/gtest.h>
#include <iostream>

typedef int64_t IntegerType;

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
    }

    // virtual void TearDown() {}
    VectorArray *v1;
};

TEST_F(ParallelGraverTest, permuteRank) {
    VectorArray *v2 = new VectorArray (*v1); //copy
    Permutation p = ParallelGraver::permute_full_rank_to_left (*v2);
    EXPECT_FALSE(*v2 == *v1);
    v2->permute(p);
    EXPECT_TRUE(*v2 == *v1);
}
