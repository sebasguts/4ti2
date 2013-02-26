#include <gtest/gtest.h>
#include <iostream>

typedef int64_t IntegerType;

#include "groebner/Feasible.h"
#include "groebner/FeasibleStream.h"
#include "groebner/LatticeBasis.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/GraverVectors.h"

using namespace _4ti2_;

class GraverVectorsTest : public ::testing::Test {
protected:
    virtual void SetUp() {
    }

    // virtual void TearDown() {}
};

TEST_F(GraverVectorsTest, create) {
    EXPECT_TRUE( 0==0 );
}
