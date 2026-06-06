#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/geometry.h>


TEST_CASE("normalized angles"){
	const std::vector<double> angles = {
		0,
		M_PI - 1e-6,
		-M_PI + 1e-6,
		M_PI/2,
		-M_PI/2,
		1.234,
		-2.345,
	};
	for(const auto &angle : angles){
		CHECK_THAT(math::normalize_angle(angle), Catch::Matchers::WithinAbs(angle, 1e-10));
		CHECK_THAT(math::normalize_angle(angle + 2*M_PI), Catch::Matchers::WithinAbs(angle, 1e-10));
		CHECK_THAT(math::normalize_angle(angle - 2*M_PI), Catch::Matchers::WithinAbs(angle, 1e-10));
		CHECK_THAT(math::normalize_angle(angle + 10 * 2*M_PI), Catch::Matchers::WithinAbs(angle, 1e-10));
		CHECK_THAT(math::normalize_angle(angle - 10 * 2*M_PI), Catch::Matchers::WithinAbs(angle, 1e-10));
	}
}


TEST_CASE("rotate vectors"){
	const math::Vector3d v0 = {0, 0, 0};
	const math::Vector3d vx = {1, 0, 0};
	const math::Vector3d vy = {0, 1, 0};
	const math::Vector3d vz = {0, 0, 1};

	SECTION("180 degrees around every axis leads back to the start"){
		const auto R = math::rotate3_xyz(M_PI, M_PI, M_PI);

		const auto r0 = (R*v0.extended()).projected();
		const auto rx = (R*vx.extended()).projected();
		const auto ry = (R*vy.extended()).projected();
		const auto rz = (R*vz.extended()).projected();
		CHECK_THAT(r0.norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((rx - vx).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((ry - vy).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((rz - vz).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
	SECTION("90 degrees extrinsic around every axis"){
		const auto R = math::rotate3_xyz(M_PI/2, M_PI/2, M_PI/2);

		const auto r0 = (R*v0.extended()).projected();
		const auto rx = (R*vx.extended()).projected();
		const auto ry = (R*vy.extended()).projected();
		const auto rz = (R*vz.extended()).projected();
		CHECK_THAT(r0.norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((rx - (-1.0*vz)).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((ry - vy).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((rz - vx).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
	SECTION("90 degrees intrinsic around every axis"){
		const auto R = math::rotate3_rpy(M_PI/2, M_PI/2, M_PI/2);

		const auto r0 = (R*v0.extended()).projected();
		const auto rx = (R*vx.extended()).projected();
		const auto ry = (R*vy.extended()).projected();
		const auto rz = (R*vz.extended()).projected();
		CHECK_THAT(r0.norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((rx - vz).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((ry - (-1.0*vy)).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((rz - vx).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
}


TEST_CASE("rotation equivalences"){
	const auto id4 = math::Matrix<double, 4, 4>::I();

	SECTION("rotation equivalences 1"){
		const math::Vector3d diag = math::Vector3d({1, 1, -1}).normalized();
		const auto R = math::rotate3_axis(diag, 2*M_PI/3);
		const auto Ri = R.inverse();
		REQUIRE_THAT((R * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));

		const auto R_1 = math::rotate3_xyz(0.0, M_PI/2, -M_PI/2);
		const auto R_2 = math::rotate3_rpy(M_PI/2, 0.0, -M_PI/2);
		const auto R_3 = math::rotate3_quaternion(math::Vector4d({1, 1, 1, -1}));
		const auto R_4 = math::rotate3_quaternion(math::Vector3d({0.5, 0.5, -0.5}));

		CHECK_THAT((R_1 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((R_2 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((R_3 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((R_4 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}

	SECTION("rotation equivalences 2"){
		const math::Vector3d diag = math::Vector3d({1, 1, 1}).normalized();
		const auto R = math::rotate3_axis(diag, -2*M_PI/3);
		const auto Ri = R.inverse();
		REQUIRE_THAT((R * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));

		const auto R_1 = math::rotate3_xyz(M_PI/2, -M_PI/2, M_PI);
		const auto R_2 = math::rotate3_rpy(-M_PI/2, 0.0, -M_PI/2);
		const auto R_3 = math::rotate3_quaternion(math::Vector4d({1, -1, -1, -1}));
		const auto R_4 = math::rotate3_quaternion(math::Vector3d({-0.5, -0.5, -0.5}));

		CHECK_THAT((R_1 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((R_2 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((R_3 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
		CHECK_THAT((R_4 * Ri - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}

	SECTION("extract axis and angle"){
		const math::Vector3d diag = math::Vector3d({21.3, 17.71, 99.99}).normalized();
		const double rot = 4*M_PI / 3;
		const auto R = math::rotate3_axis(diag, rot);
		const auto R_i = R.inverse();
		REQUIRE_THAT((R * R_i - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));

		double rot_2;
		const math::Vector3d diag_2 = math::get_axis(R).normalized(rot_2);

		const bool same = Catch::Matchers::WithinAbs(0, 1e-10).match((diag - diag_2).norm()) && Catch::Matchers::WithinAbs(0, 1e-10).match(math::angle_diff(rot, rot_2));
		const bool inv = Catch::Matchers::WithinAbs(0, 1e-10).match((diag + diag_2).norm()) && Catch::Matchers::WithinAbs(0, 1e-10).match(math::normalize_angle(rot + rot_2));
		CHECK((same || inv));

		const auto R_1 = math::rotate3_axis(diag_2, rot_2);
		CHECK_THAT((R_1 * R_i - id4).norm(), Catch::Matchers::WithinAbs(0, 1e-10));
	}
}


