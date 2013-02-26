#include <gtest/gtest.h>
#include <iostream>

typedef int64_t IntegerType;

#include "groebner/Vector.h"
#include "groebner/VectorArray.h"
#include "groebner/GraverVector.h"

using namespace _4ti2_;

class GraverVectorTest : public ::testing::Test {
protected:
    virtual void SetUp() {
	va = new VectorArray (4,5,1); // all ones matrix;
	(*va)[0][0] = 1;
	(*va)[0][1] = 1;
	(*va)[0][2] = 0;
	(*va)[0][3] = 0;
	(*va)[0][4] = 0;

	(*va)[1][0] = -1;
	(*va)[1][1] = 1;
	(*va)[1][2] = 0;
	(*va)[1][3] = 0;
	(*va)[1][4] = 0;

	(*va)[2][0] = 0;
	(*va)[2][1] = 0;
	(*va)[2][2] = 17;
	(*va)[2][3] = 19;
	(*va)[2][4] = 0;

	(*va)[3][0] = 0;
	(*va)[3][1] = 0;
	(*va)[3][2] = -3;
	(*va)[3][3] = 0;
	(*va)[3][4] = 1;	    

        v1 = new GraverVector ( & (*va)[0] );
        v2 = new GraverVector ( & (*va)[1] );
	v3 = new GraverVector ( & (*va)[2] );
	v4 = new GraverVector ( & (*va)[3] );
	v5 = new GraverVector ( *v4 );
    }

    // virtual void TearDown() {}
    GraverVector *v1,*v2,*v3, *v4;
    GraverVector *v5;
    VectorArray *va;
};

TEST_F(GraverVectorTest, create) {
    EXPECT_TRUE (v4->v == v5->v); // Compare adresses
    EXPECT_TRUE ( *(v4->v) == *(v5->v) ); // Run vector comparison
}

TEST_F(GraverVectorTest, assign) {
    GraverVector w = *v2;
    EXPECT_EQ(w.v ,v2->v);
}

TEST_F(GraverVectorTest, lastEntry) {
    EXPECT_EQ(v1->last_entry(), 0 );
    EXPECT_EQ(v4->last_entry(), 1 );
}

TEST_F(GraverVectorTest, signConsistent) {
    EXPECT_TRUE (v4->is_sign_consistent(*v5));
    EXPECT_FALSE (v4->is_sign_consistent(*v3));
    EXPECT_FALSE (v1->is_sign_consistent(*v2));
    EXPECT_TRUE (v1->is_sign_consistent(*v3));

    // Bug from the code: 0 -1  2 -1  + -4  0  0  1;
    // Those two vectors are sign consistent since we only look at the first n-1 signs!
    Vector v(4,0),w(4,0);
    v[0]=0; v[1]=-1; v[2]=2; v[3]=-1;
    w[0]=-4; w[1]=0; w[2]=0; w[3]=1;
    GraverVector gv(&v);
    GraverVector gw(&w);
    EXPECT_TRUE(gv.is_applicable(gw));
}
