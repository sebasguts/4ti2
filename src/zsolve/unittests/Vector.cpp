#include <gtest/gtest.h>

#include "../Vector.hpp"

using namespace _4ti2_zsolve_;

TEST(Vector, consistency) {
    size_t *v = create_vector <size_t> (17);
    EXPECT_TRUE(check_vector_consistency <size_t> (v, 17));
}

TEST(Vector, values) {
    size_t *v = create_vector <size_t> (17,3);
    EXPECT_EQ(3,v[4]);
}
