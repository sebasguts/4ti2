#include <gtest/gtest.h>

#include "../VectorArray.hpp"

typedef std::vector<int> T;

class VectorArrayTest : public ::testing::Test {
protected:
    virtual void SetUp() {
	t0 = new _4ti2_zsolve_::VectorArray<T> (0);
	t1 = new _4ti2_zsolve_::VectorArray<T> (10);
	t2 = new _4ti2_zsolve_::VectorArray<T> (5,10);
    }

    // virtual void TearDown() {}
    
    _4ti2_zsolve_::VectorArray<T>* t0;
    _4ti2_zsolve_::VectorArray<T>* t1;
    _4ti2_zsolve_::VectorArray<T>* t2;
};

TEST_F(VectorArrayTest, sizes) {
    EXPECT_EQ(0, t0->height());
    EXPECT_EQ(0, t0->width());
    EXPECT_EQ(5, t2->height());
    EXPECT_EQ(10, t2->width());
}

TEST_F(VectorArrayTest, sizeAfterInsert) {
    T v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    t2->append_vector(&v);
    EXPECT_EQ(6, t2->height());
}
