/*
 * AllTests.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: bryan
 */
#include <gtest/gtest.h>
#include "VectorTests.cpp"
#include "DiskTests.cpp"

int main (int argc, char **argv){
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
