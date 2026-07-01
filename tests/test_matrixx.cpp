#include <catch2/catch_test_macros.hpp>

#include <math/Matrix.h>
#include <math/MatrixX.h>


TEST_CASE("MatrixX to Matrix"){
	math::MatrixX<int> mx(4, 5);
	mx = {
		1, 2, 3, 4, 5,
		1, 2, 3, 4, 5,
		1, 2, 3, 4, 5,
		1, 2, 3, 4, 5,
	};
	const math::Matrix<int, 4, 5> mnx = mx.get_data();
	for(size_t i=0; i<4; i++){
		for(size_t j=0; j<5; j++){
			CHECK(mx(i, j) == mnx(i, j));
		}
	}
}


TEST_CASE("MatrixX Dimension Mismatches and Exceptions") {
	math::MatrixX<double> m1(2, 3);
	math::MatrixX<double> m2(3, 4);
	math::MatrixX<double> m3(2, 2);

	SECTION("Invalid Initialization") {
		CHECK_THROWS_AS((m1 = {1.0, 2.0, 3.0}), std::logic_error);
	}

	SECTION("Out of Bounds Indexing") {
		CHECK_THROWS_AS(m1(2, 0), std::logic_error);
		CHECK_THROWS_AS(m1(0, 3), std::logic_error);
	}

	SECTION("Invalid Matrix Multiplication Dimensions") {
		CHECK_THROWS_AS(m3 * m1, std::logic_error);
	}

	SECTION("Invalid Addition Dimensions") {
		CHECK_THROWS_AS(m1 += m2, std::logic_error);
	}
}

TEST_CASE("MatrixX Operations") {
	math::MatrixX<int> m1(2, 3);
	m1 = {
		1, 2, 3,
		4, 5, 6
	};

	SECTION("Transposition") {
		auto mT = m1.transposed();
		REQUIRE(mT.rows() == 3);
		REQUIRE(mT.cols() == 2);
		CHECK(mT(0, 0) == 1);
		CHECK(mT(0, 1) == 4);
		CHECK(mT(1, 0) == 2);
		CHECK(mT(1, 1) == 5);
		CHECK(mT(2, 0) == 3);
		CHECK(mT(2, 1) == 6);
	}

	SECTION("Scalar Multiplication and Division") {
		auto m2 = m1 * 2;
		CHECK(m2(0, 0) == 2);
		CHECK(m2(1, 2) == 12);

		auto m3 = m2 / 2;
		CHECK(m3 == m1);
	}
}


