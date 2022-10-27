#define CATCH_CONFIG_MAIN

#include <catch/catch.hpp>

#include <iostream>

#include "mpreal.h"

// Test omega function.

TEST_CASE("An example test", "example1") {

    int i = 1;
    int j = 1;

    std::cout << "This is an example test" << std::endl;

    REQUIRE(i == j);

}
