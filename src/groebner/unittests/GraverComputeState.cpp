#include <gtest/gtest.h>
#include <iostream>

typedef int64_t IntegerType;

#include "groebner/Feasible.h"
#include "groebner/FeasibleStream.h"
#include "groebner/LatticeBasis.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/GraverComputeState.h"

using namespace _4ti2_;

class GraverComputeStateTest : public ::testing::Test {
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

TEST_F(GraverComputeStateTest, create) {
    GraverComputeState s (*v2);
    // s.MakeGraverBasisWithZSolve();
}
