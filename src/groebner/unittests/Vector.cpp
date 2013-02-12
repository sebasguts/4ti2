#include <gtest/gtest.h>

typedef int64_t IntegerType;

#include "groebner/Vector.h"

using namespace _4ti2_;

class VectorTest : public ::testing::Test {
protected:
    virtual void SetUp() {
	v1 = new Vector;
	v2 = new Vector (10);
	v3 = new Vector (10, 17);
	v4 = new Vector (*v3);
	v5 = new Vector (-(*v4));
    }

    // virtual void TearDown() {}
    Vector *v1;
    Vector *v2;
    Vector *v3;
    Vector *v4;
    Vector *v5;
};

TEST_F(VectorTest, Construction) {
    EXPECT_TRUE(*v3 == *v4);
}

TEST_F(VectorTest, Norm) {
    EXPECT_EQ(v4->norm(), 170);
    EXPECT_EQ(v4->norm(), v5->norm());
    EXPECT_EQ(v4->norm(1), 17);
}

TEST(Vector, move) {
    Vector *v = new Vector(10,3);
    Vector w = std::move(*v); // no copy called
    delete v; // <- still no exception;
    Vector x = Vector (10,4); // Only single construction, no copy
}
