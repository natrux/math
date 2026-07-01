#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/Matrix.h>
#include <math/Vector.h>


template<class T, size_t N>
math::Vector<T, N> solve_equations(const math::Matrix<T, N, N> &lhs, const math::Vector<T, N> &rhs){
	math::Matrix<T, N, N+1> H;
	H.set(0, 0, lhs);
	H.set(0, N, rhs);
	H.gauss_jordan();
	for(size_t i=0; i<3; i++){
		REQUIRE(H(i, i) != 0);
		H.set(i, 0, H.get_row(i) / H(i, i));
	}
	const auto sol = H.get_column(3);
	const auto R = rhs - lhs * sol;
	CHECK_THAT(R.norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	return sol;
}


TEST_CASE("determinants and inversions"){
	const auto id0 = math::Matrix<double, 0, 0>::I();
	const auto id1 = math::Matrix<double, 1, 1>::I();
	const auto id2 = math::Matrix<double, 2, 2>::I();
	const auto id3 = math::Matrix<double, 3, 3>::I();
	const auto id4 = math::Matrix<double, 4, 4>::I();
	const auto id5 = math::Matrix<double, 5, 5>::I();
	const auto id6 = math::Matrix<double, 6, 6>::I();

	SECTION("0x0"){
		const math::Matrix<double, 0, 0> m0;
		const math::Matrix<int, 0, 0> m0_int = m0;
		const auto m0_1 = m0.inverse();
		CHECK(m0.determinant() == 1);
		CHECK(m0_int.determinant() == 1);
		CHECK(m0 * m0_1 == id0);
	}
	SECTION("1x1"){
		const math::Matrix<double, 1, 1> m1 = {
			5,
		};
		const math::Matrix<int, 1, 1> m1_int = m1;
		const auto m1_1 = m1.inverse();
		CHECK(m1.determinant() == 5);
		CHECK(m1_int.determinant() == 5);
		CHECK_THAT((m1 * m1_1 - id1).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
	SECTION("2x2"){
		const math::Matrix<double, 2, 2> m2 = {
			1, 2,
			3, 4,
		};
		const math::Matrix<int, 2, 2> m2_int = m2;
		const auto m2_1 = m2.inverse();
		CHECK_THAT(m2.determinant(), Catch::Matchers::WithinAbs(-2, 1e-10));
		CHECK(m2_int.determinant() == m2.determinant());
		CHECK_THAT((m2 * m2_1 - id2).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
	SECTION("3x3"){
		const math::Matrix<double, 3, 3> m3 = {
			2, -1, 0,
			-1, 2, -1,
			0, -1, 2,
		};
		const math::Matrix<int, 3, 3> m3_int = m3;
		const auto m3_1 = m3.inverse();
		CHECK_THAT(m3.determinant(), Catch::Matchers::WithinAbs(4, 1e-10));
		CHECK(m3_int.determinant() == 4);
		CHECK_THAT((m3 * m3_1 - id3).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
	SECTION("4x4"){
		const math::Matrix<double, 4, 4> m4 = {
			1, 1, 1, -1,
			1, 1, -1, 1,
			1, -1, 1, 1,
			-1, 1, 1, 1,
		};
		const math::Matrix<int, 4, 4> m4_int = m4;
		const auto m4_1 = m4.inverse();
		CHECK_THAT(m4.determinant(), Catch::Matchers::WithinAbs(-16, 1e-10));
		CHECK(m4_int.determinant() == -16);
		CHECK_THAT((m4 * m4_1 - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
	SECTION("5x5"){
		const math::Matrix<double, 5, 5> m5 = {
			4, -1, -4, -3,  4,
			2, -1,  2, -2, -1,
			3, -2, -1, -3,  1,
			2,  1,  1, -1,  1,
			1,  0,  2, -1, -1,
		};
		const math::Matrix<int, 5, 5> m5_int = m5;
		const auto m5_1 = m5.inverse();
		CHECK_THAT(m5.determinant(), Catch::Matchers::WithinAbs(1, 1e-10));
		CHECK(m5_int.determinant() == 1);
		CHECK_THAT((m5 * m5_1 - id5).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
	SECTION("6x6"){
		const math::Matrix<double, 6, 6> m6 = {
			1, 1, 1, 1, 1, -1,
			1, 1, 1, 1, -1, 1,
			1, 1, 1, -1, 1, 1,
			1, 1, -1, 1, 1, 1,
			1, -1, 1, 1, 1, 1,
			-1, 1, 1, 1, 1, 1,
		};
		const math::Matrix<int, 6, 6> m6_int = m6;
		const auto m6_1 = m6.inverse();
		CHECK_THAT(m6.determinant(), Catch::Matchers::WithinAbs(128, 1e-10));
		CHECK(m6_int.determinant() == 128);
		CHECK_THAT((m6 * m6_1 - id6).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
}


TEST_CASE("equation systems"){
	solve_equations<double, 3>(
		{
			1, 1, 1,
			4, 2, 1,
			9, 3, 1,
		},
		{0, 1, 3}
	);
	solve_equations<double, 3>(
		{
			7, 3, 2,
			-5, 2, 1,
			0, 12, 1,
		},
		{-2, 9, 11}
	);
}


TEST_CASE("Matrix Error States & Integer Gauss-Jordan") {
	SECTION("Singular Matrix Inverse Detection") {
		math::Matrix2d singular = {
			1.0, 2.0,
			2.0, 4.0
		};
		CHECK_FALSE(singular.is_invertible());
		CHECK_THROWS_AS(singular.inverse(), std::logic_error);
	}

	SECTION("Index Bounds Validation") {
		math::Matrix2d valid;
		CHECK_THROWS_AS(valid(2, 0), std::logic_error);
		CHECK_THROWS_AS(valid[4], std::logic_error);
	}

	SECTION("Integer Matrix Elimination Precision") {
		math::Matrix<int, 2, 2> int_matrix = {
			3, 4,
			2, 3
		};
		int_matrix.gauss_jordan();

		CHECK(int_matrix(0, 0) != 0);
		CHECK(int_matrix(0, 1) == 0);
		CHECK(int_matrix(1, 0) == 0);
		CHECK(int_matrix(1, 1) != 0);
	}
}

