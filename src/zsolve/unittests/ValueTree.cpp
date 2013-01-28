#include <gtest/gtest.h>

bool vt () { return true; };

TEST(ValueTree, create) {
    EXPECT_EQ(true, vt());
    EXPECT_FALSE(!(vt()));
}

TEST(ValueTree, destroy) {
    EXPECT_EQ(true, vt());
    EXPECT_FALSE(!(vt()));
}

