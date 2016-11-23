#include <gtest/gtest.h>
#include "Vector.h"

class VectorTest : public ::testing::Test{
	protected:
		Vector v0;
		Vector v1 = Vector(1.0,1.0,1.0);
		Vector v2 = Vector(2.0,2.0,2.0);
		Vector v3 = Vector(3.0,3.0,3.0);
};

TEST_F(VectorTest,VectorImplicitOK){
	ASSERT_EQ(v0.x(),0.0);
	ASSERT_EQ(v0.y(),0.0);
	ASSERT_EQ(v0.phi(),0.0);
}
TEST_F(VectorTest,VectorExplicitOK){
	ASSERT_EQ(v1.x(),1.0);
	ASSERT_EQ(v1.y(),1.0);
	ASSERT_EQ(v1.phi(),1.0);
	ASSERT_EQ(v2.x(),2.0);
	ASSERT_EQ(v2.y(),2.0);
	ASSERT_EQ(v2.phi(),2.0);
	ASSERT_EQ(v3.x(),3.0);
	ASSERT_EQ(v3.y(),3.0);
	ASSERT_EQ(v3.phi(),3.0);
}
TEST_F(VectorTest,VectorAdditionOK){
	Vector v = v1;
	v+=v2;
	ASSERT_EQ(v.x(),v3.x());
	ASSERT_EQ(v.x(),v3.x());
	ASSERT_EQ(v.x(),v3.x());
}
TEST_F(VectorTest,VectorSubtractionOK){
	Vector v = v3;
	v-=v2;
	ASSERT_EQ(v.x(),v1.x());
	ASSERT_EQ(v.x(),v1.x());
	ASSERT_EQ(v.x(),v1.x());
}
TEST_F(VectorTest,MultiplicationByScalarOK){
	Vector v = v1;
	v*=3.0;
	ASSERT_EQ(v.x(),v3.x());
	ASSERT_EQ(v.x(),v3.x());
	ASSERT_EQ(v.x(),v3.x());
}

