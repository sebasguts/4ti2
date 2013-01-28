#include <gtest/gtest.h>

bool testfcn() { return true; }

TEST(Nothing, ideal) {
    EXPECT_EQ(true, testfcn());
    EXPECT_FALSE(!(testfcn()));
}

