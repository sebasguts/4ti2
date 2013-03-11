#include <gtest/gtest.h>
#include <iostream>

typedef int64_t IntegerType;

#include <string>

#include "groebner/VectorArrayStream.h"
#include "groebner/LatticeBasis.h"

using namespace _4ti2_;

class LatticeBasisTest : public ::testing::Test {
protected:
    virtual void SetUp() {
    }

    // virtual void TearDown() {}
    VectorArray A, R;
};


TEST_F(LatticeBasisTest, IdentitySolve) {
    A = VectorArray (3,3);
    std::stringstream As;
    As << "1 0 0 0 1 0 0 0 1";
    R = VectorArray (2,3);
    std::stringstream Rs;
    Rs << "1 2 3 4 5 6";
    As >> A;
    Rs >> R;

    VectorArray S (2,3);
    solve (A, R, S);
    EXPECT_TRUE (S == R);
}

TEST_F(LatticeBasisTest, 44solve) {
    A = VectorArray (4,4);
    std::stringstream As;
    As << "1 2 3 4 5 6 7 8 17 3 4 5 6 7 8 11";
    As >> A;
    R = VectorArray (2,4);
    std::stringstream Rs;
    Rs << "11 35 32 37 -70 -414 -1671 -494";
    Rs >> R;
    VectorArray S (2,4);
    std::stringstream Ss;
    Ss << "1 3 4 -2 -101 17 -5 3";
    Ss >> S;
    VectorArray S2 (2,4);
    solve (A, R, S2);
    EXPECT_TRUE (S == S2);
}

