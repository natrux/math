#include <math/math.h>

#include <iostream>


template<class T, size_t N, size_t M>
void print_matrix(const math::Matrix<T, N, M> &matrix){
	for(size_t i=0; i<N; i++){
		for(size_t j=0; j<M; j++){
			if(j > 0){
				std::cout << " ";
			}
			std::cout << matrix(i, j);
		}
		std::cout << std::endl;
	}
}


template<class T, size_t N>
math::Vector<T, N> solve_equations(const math::Matrix<T, N, N> &lhs, const math::Vector<T, N> &rhs){
	math::Matrix<T, N, N+1> H;
	H.set(0, 0, lhs);
	H.set(0, N, rhs);
	H.gauss_jordan();
	for(size_t i=0; i<3; i++){
		H.set(i, 0, H.get_row(i) / H(i, i));
	}
	const auto sol = H.get_column(3);
	const auto R = rhs - lhs * sol;
	if(R.norm() >= 1e-10){
		throw std::runtime_error("Equation system wrong solution");
	}
	return sol;
}


static void test(){
	math::Matrix<float, 5, 5> five;
	math::Matrix<double, 4, 4> four;
	math::Matrix<double, 3, 3> three;
	math::Matrix<int, 2, 2> two;
	math::Matrix<double, 1, 1> one;
	math::Matrix<double, 0, 0> zero;

	const auto id0 = math::Matrix<double, 0, 0>::I();
	const auto id1 = math::Matrix<double, 1, 1>::I();
	const auto id2 = math::Matrix<double, 2, 2>::I();
	const auto id3 = math::Matrix<double, 3, 3>::I();
	const auto id4 = math::Matrix<double, 4, 4>::I();
	const auto id5 = math::Matrix<double, 5, 5>::I();
	const auto id6 = math::Matrix<double, 6, 6>::I();

	const math::Matrix<double, 0, 0> m0;
	const auto m0_1 = m0.inverse();
	if(m0 * m0_1 != id0){
		throw std::runtime_error("failed to invert m0");
	}

	const math::Matrix<double, 1, 1> m1 = {
		5,
	};
	if(m1.determinant() != 5){
		throw std::runtime_error("wrong determinant for m1");
	}
	const auto m1_1 = m1.inverse();
	if((m1 * m1_1 - id1).norm() > 1e-10){
		throw std::runtime_error("failed to invert m1");
	}

	const math::Matrix<double, 2, 2> m2 = {
		1, 2,
		3, 4,
	};
	if(m2.determinant() != -2){
		throw std::runtime_error("wrong determinant for m2");
	}
	const auto m2_1 = m2.inverse();
	if((m2 * m2_1 - id2).norm() > 1e-10){
		throw std::runtime_error("failed to invert m2");
	}

	const math::Matrix<double, 3, 3> m3 = {
		2, -1, 0,
		-1, 2, -1,
		0, -1, 2,
	};
	if(m3.determinant() != 4){
		throw std::runtime_error("wrong determinant for m3");
	}
	const auto m3_1 = m3.inverse();
	if((m3 * m3_1 - id3).norm() > 1e-10){
		throw std::runtime_error("failed to invert m3");
	}

	const math::Matrix<double, 4, 4> m4 = {
		1, 1, 1, -1,
		1, 1, -1, 1,
		1, -1, 1, 1,
		-1, 1, 1, 1,
	};
	if(m4.determinant() != -16){
		throw std::runtime_error("wrong determinant for m4");
	}
	const auto m4_1 = m4.inverse();
	if((m4 * m4_1 - id4).norm() > 1e-10){
		throw std::runtime_error("failed to invert m4");
	}

	const math::Matrix<double, 5, 5> m5 = {
		4, -1, -4, -3,  4,
		2, -1,  2, -2, -1,
		3, -2, -1, -3,  1,
		2,  1,  1, -1,  1,
		1,  0,  2, -1, -1,
	};
	if(m5.determinant() != 1){
		throw std::runtime_error("wrong determinant for m5");
	}
	{
		const math::Matrix<int, 5, 5> m5_int = m5;
		if(m5_int.determinant() != m5.determinant()){
			throw std::runtime_error("determinant differs between types");
		}
	}
	const auto m5_1 = m5.inverse();
	if((m5 * m5_1 - id5).norm() > 1e-10){
		throw std::runtime_error("failed to invert m5");
	}

	const math::Matrix<double, 6, 6> m6 = {
		1, 1, 1, 1, 1, -1,
		1, 1, 1, 1, -1, 1,
		1, 1, 1, -1, 1, 1,
		1, 1, -1, 1, 1, 1,
		1, -1, 1, 1, 1, 1,
		-1, 1, 1, 1, 1, 1,
	};
	if(m6.determinant() != 128){
		throw std::runtime_error("wrong determinant for m6");
	}
	const auto m6_1 = m6.inverse();
	if((m6 * m6_1 - id6).norm() > 1e-10){
		throw std::runtime_error("failed to invert m6");
	}

	math::Vector3d vector = {1, 0, 0};
	vector.x();
	vector = {0, 1, 0};
	const math::Vector3d vector_2 = vector;
	vector_2.z();
	const auto vector_3 = vector.dotted(vector_2);
	vector_3[1];
	const math::Vector2d vector2 = {-2, -2};
	if(std::abs(vector2.atan2() - (-3*M_PI/4))  > 1e-10){
		throw std::runtime_error("Wrong angle: " + std::to_string(vector2.atan2()));
	}

	(math::rotate3_rpy(1.0, 2.0, 3.0) * vector_2.extended()).projected();

	solve_equations<double, 3>(
		{
			1, 1, 1,
			4, 2, 1,
			9, 3, 1,
		},
		{0, 1, 3}
	);

	math::Matrix<double, 2, 2> uninitialized = {};
	std::array<double, 4> bla;
	uninitialized = bla;

	const math::Vector3i cross_1 = {3, 5, 11};
	const math::Vector3i cross_2 = {23, 42, -99};
	if(cross_1.crossed(cross_2) != cross_1.crossed() * cross_2){
		throw std::runtime_error("cross product matrix is wrong");
	}

	const math::Matrix<int, 3, 4> crosses_1 = {
		12, 22, 32, 43,
		88, 99, 71, 9,
		-4, -77, -21, -1234,
	};
	const math::Matrix<int, 3, 4> crosses_2 = {
		-99, 99, -4, 44,
		33, 333, 3333, 0,
		17, 71, 23, 42,
	};
	const auto crossesed = crosses_1.crossed(crosses_2);
	for(size_t i=0; i<4; i++){
		if(crosses_1.get_column(i).crossed(crosses_2.get_column(i)) != crossesed.get_column(i)){
			throw std::runtime_error("Matrix cross product extension wrong");
		}
	}

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
			if(mx(i, j) != mnx(i, j)){
				throw std::runtime_error("Converting MatrixX to Matrix failed");
			}
		}
	}

	{
		const math::Vector3d diag = math::Vector3d({1, 1, -1}).normalized();
		const auto R = math::rotate3_axis(diag, 2*M_PI/3);
		const auto Ri = R.inverse();

		const auto R_1 = math::rotate3_xyz(0.0, M_PI/2, -M_PI/2);
		const auto R_2 = math::rotate3_rpy(M_PI/2, 0.0, -M_PI/2);
		const auto R_3 = math::rotate3_quaternion(math::Vector4d({1, 1, 1, -1}));
		const auto R_4 = math::rotate3_quaternion(math::Vector3d({0.5, 0.5, -0.5}));

		if(
			(R * Ri - id4).norm() > 1e-10 ||
			(R_1 * Ri - id4).norm() > 1e-10 ||
			(R_2 * Ri - id4).norm() > 1e-10 ||
			(R_3 * Ri - id4).norm() > 1e-10 ||
			(R_4 * Ri - id4).norm() > 1e-10 ||
			false
		){
			throw std::runtime_error("Rotations not equivalent");
		}
	}
	{
		const math::Vector3d diag = math::Vector3d({1, 1, 1}).normalized();
		const auto R = math::rotate3_axis(diag, -2*M_PI/3);
		const auto Ri = R.inverse();

		const auto R_1 = math::rotate3_xyz(M_PI/2, -M_PI/2, M_PI);
		const auto R_2 = math::rotate3_rpy(-M_PI/2, 0.0, -M_PI/2);
		const auto R_3 = math::rotate3_quaternion(math::Vector4d({1, -1, -1, -1}));
		const auto R_4 = math::rotate3_quaternion(math::Vector3d({-0.5, -0.5, -0.5}));

		if(
			(R * Ri - id4).norm() > 1e-10 ||
			(R_1 * Ri - id4).norm() > 1e-10 ||
			(R_2 * Ri - id4).norm() > 1e-10 ||
			(R_3 * Ri - id4).norm() > 1e-10 ||
			(R_4 * Ri - id4).norm() > 1e-10 ||
			false
		){
			throw std::runtime_error("Rotations not equivalent");
		}
	}
	{
		const math::Vector3d diag = math::Vector3d({21.3, 17.71, 99.99}).normalized();
		const double rot = 4*M_PI / 3;
		const auto R = math::rotate3_axis(diag, rot);
		const auto R_i = R.inverse();

		double rot_2;
		const math::Vector3d diag_2 = math::get_axis(R).normalized(rot_2);

		if(
			!((diag - diag_2).norm() < 1e-10 && std::abs(math::angle_diff(rot, rot_2)) < 1e-10) &&
			!((diag + diag_2).norm() < 1e-10 && std::abs(math::normalize_angle(rot + rot_2)) < 1e-10)
		){
			throw std::runtime_error("getting rotation axis/angle failed");
		}

		const auto R_1 = math::rotate3_axis(diag_2, rot_2);
		if((R_1 * R_i - id4).norm() > 1e-10){
			throw std::runtime_error("using calculated rotation axis/angle failed");
		}
	}
	if(
		math::gcd(3528, 3780) != 252 ||
		math::gcd(33, 44, 55, 66, 77, 88) != 11 ||
		math::gcd(-10, -15) != 5 ||
		false
	){
		throw std::runtime_error("wrong gcd");
	}
	if(
		math::lcm(3528, 3780) != 52920 ||
		math::lcm(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) != 2520 ||
		math::lcm(-10, -15) != 30 ||
		false
	){
		throw std::runtime_error("wrong lcm");
	}
}


int main(int /*argc*/, char **/*argv*/){
	try{
		test();
	}catch(const std::exception &err){
		std::cerr << "Error: " << err.what() << std::endl;
		return 1;
	}
	return 0;
}
