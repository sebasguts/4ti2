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

TEST_F(VectorArrayTest, consistency) {
    EXPECT_FALSE(t0->check_consistency());
    EXPECT_TRUE(t1->check_consistency());
    EXPECT_TRUE(t2->check_consistency());
}

TEST_F(VectorArrayTest, sizes) {
    EXPECT_EQ(0, t0->height());
    EXPECT_EQ(0, t0->width());
    EXPECT_EQ(5, t2->height());
    EXPECT_EQ(10, t2->width());
    EXPECT_FALSE(t0->check_consistency());
    EXPECT_TRUE(t1->check_consistency());
    EXPECT_TRUE(t2->check_consistency());
}

TEST_F(VectorArrayTest, sizeAfterModify) {
    T* v = new T({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    t2->append_vector(v);
    EXPECT_EQ(6, t2->height());
    t2->remove_unsorted(1);
    EXPECT_EQ(5, t2->height());
    EXPECT_TRUE(t2->check_consistency());
    delete v;
}

TEST_F(VectorArrayTest, bracket) {
    T v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    t2->append_vector(&v);
    EXPECT_EQ( (*(*t2)[5]) [3] , 4);
}

TEST_F(VectorArrayTest, append_returns_last_valid_index) {
    T* v = new T({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    EXPECT_EQ(5, t2->append_vector(v));
    EXPECT_EQ(6, t2->num_vectors());
}
