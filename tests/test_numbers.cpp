#include <catch2/catch_test_macros.hpp>

#include <math/numbers.h>


TEST_CASE("greatest common divisor"){
	CHECK(math::gcd(3528, 3780) == 252);
	CHECK(math::gcd(33, 44, 55, 66, 77, 88) == 11);
	CHECK(math::gcd(-10, -15) == 5);
}


TEST_CASE("least common multiple"){
	CHECK(math::lcm(3528, 3780) == 52920);
	CHECK(math::lcm(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) == 2520);
	CHECK(math::lcm(-10, -15) == 30);
}


