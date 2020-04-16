#include "gtest/gtest.h"
  
#include "sophia_interface.h"

TEST(First, debug) {
	auto sopi = sophia_interface();

	EXPECT_THROW(sopi.debug("some string", true), std::runtime_error);
}

int main(int argc, char **argv) {
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
}
