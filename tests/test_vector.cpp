#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/Vector.h>


TEST_CASE("vector angle"){
	const math::Vector2d vector2 = {-2, -2};
	CHECK_THAT(vector2.atan2(), Catch::Matchers::WithinAbs(-3*M_PI/4, 1e-10));
}


TEST_CASE("cross product"){
	const math::Vector3i v1 = {4, -3, 9};
	const math::Vector3i v2 = {7, 12, -1};
	const math::Vector3i result = {-105, 67, 69};
	CHECK(v1.crossed(v2) == result);
}


TEST_CASE("cross product matrix"){
	const math::Vector3i cross_1 = {3, 5, 11};
	const math::Vector3i cross_2 = {23, 42, -99};
	CHECK(cross_1.crossed(cross_2) == cross_1.crossed() * cross_2);
}


TEST_CASE("matrix cross product matrix"){
	const math::Matrix<int, 3, 4> m_1 = {
		12, 22, 32, 43,
		88, 99, 71, 9,
		-4, -77, -21, -1234,
	};
	const math::Matrix<int, 3, 4> m_2 = {
		-99, 99, -4, 44,
		33, 333, 3333, 0,
		17, 71, 23, 42,
	};
	const auto cross = m_1.crossed(m_2);
	for(size_t i=0; i<4; i++){
		CHECK(m_1.get_column(i).crossed(m_2.get_column(i)) == cross.get_column(i));
	}
}


