#include <gtest/gtest.h>
#include "disk.h"

class DiskTest : public ::testing::Test{
	protected:
		disk d0;
		disk d1;
		disk d2;
		disk d3;
};

TEST_F(DiskTest,DiskImplicitOK){
	ASSERT_EQ(d0.x(),0.0);
	ASSERT_EQ(d0.y(),0.0);
	ASSERT_EQ(d0.phi(),0.0);
}

TEST_F(DiskTest,DiskExplicitOK){
	d0.pos() = 1.0,1.0,1.0;
	ASSERT_EQ(d0.x(),1.0);
	ASSERT_EQ(d0.y(),1.0);
	ASSERT_EQ(d0.phi(),1.0);
}





