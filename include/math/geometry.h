#pragma once

#include <math/Matrix.h>
#include <math/Vector.h>


namespace math{


inline double normalize_angle(double angle){
	auto m = std::fmod(angle, 2*M_PI);
	if(m > M_PI){
		m -= 2*M_PI;
	}else if(m < -M_PI){
		m += 2*M_PI;
	}
	return m;
}


inline double angle_diff(double a, double b){
	return normalize_angle(a - b);
}


// 2D

template<class T>
Matrix3<T> translate2(T x, T y){
	return {
		1, 0, x,
		0, 1, y,
		0, 0, 1,
	};
}

template<class T>
Matrix3<T> translate2(const Vector2<T> &vec){
	return translate2(vec.x(), vec.y());
}

template<class T>
Matrix3d rotate2(T r){
	return {
		std::cos(r), -std::sin(r), 0,
		std::sin(r),  std::cos(r), 0,
		0,            0,           1,
	};
}

template<class T>
Matrix2d rotation2(T r){
	return {
		std::cos(r), -std::sin(r),
		std::sin(r),  std::cos(r),
	};
}

template<class T>
Matrix3d transform2(T x, T y, T r){
	return translate2(x, y) * rotate2(r);
}

template<class T>
Matrix3d transform2(const Vector2<T> &trans, T r){
	return transform2(trans.x(), trans.y(), r);
}


// 2.5D

template<class T>
Matrix4<T> translate25(T x, T y){
	return {
		1, 0, 0, x,
		0, 1, 0, y,
		0, 0, 1, 0,
		0, 0, 0, 1,
	};
}

template<class T>
Matrix4d rotate25(T r){
	return {
		std::cos(r), -std::sin(r), 0, 0,
		std::sin(r),  std::cos(r), 0, 0,
		0,            0,           1, r,
		0,            0,           0, 1,
	};
}

template<class T>
Matrix4d transform25(T x, T y, T r){
	return translate25(x, y) * rotate25(r);
}

template<class T>
Matrix4d transform25(const Vector3<T> &vec){
	return transform25(vec.x(), vec.y(), vec.z());
}


// 3D

template<class T>
Matrix4<T> translate3(T x, T y, T z){
	return {
		1, 0, 0, x,
		0, 1, 0, y,
		0, 0, 1, z,
		0, 0, 0, 1,
	};
}

template<class T>
Matrix4<T> translate3(const Vector3<T> &vec){
	return translate3(vec.x(), vec.y(), vec.z());
}


template<class T>
Matrix4d rotate3_x(T angle){
	return {
		1, 0,                0,               0,
		0, std::cos(angle), -std::sin(angle), 0,
		0, std::sin(angle),  std::cos(angle), 0,
		0, 0,                0,               1,
	};
}

template<class T>
Matrix4d rotate3_y(T angle){
	return {
		 std::cos(angle), 0, std::sin(angle), 0,
		 0,               1, 0,               0,
		-std::sin(angle), 0, std::cos(angle), 0,
		 0,               0, 0,               1,
	};
}

template<class T>
Matrix4d rotate3_z(T angle){
	return {
		std::cos(angle), -std::sin(angle), 0, 0,
		std::sin(angle),  std::cos(angle), 0, 0,
		0,                0,               1, 0,
		0,                0,               0, 1,
	};
}

// extrinsic rotations around x, y and z axis
template<class T>
Matrix4d rotate3_xyz(T rx, T ry, T rz){
	return rotate3_z(rz) * rotate3_y(ry) * rotate3_x(rx);
}

// extrinsic rotations around z, y and x axis
template<class T>
Matrix4d rotate3_zyx(T rx, T ry, T rz){
	return rotate3_x(rx) * rotate3_y(ry) * rotate3_z(rz);
}

// intrinsic rotations around roll, pitch and yaw axis
template<class T>
Matrix4d rotate3_rpy(T roll, T pitch, T yaw){
	return rotate3_zyx(roll, pitch, yaw);
}

// intrinsic rotations around yaw, pitch and roll axis
template<class T>
Matrix4d rotate3_ypr(T roll, T pitch, T yaw){
	return rotate3_xyz(roll, pitch, yaw);
}

template<class T>
Matrix4d rotate3_xyz(const Vector3<T> &vec){
	return rotate3_xyz(vec.x(), vec.y(), vec.z());
}

template<class T>
Matrix4d rotate3_zyx(const Vector3<T> &vec){
	return rotate3_zyx(vec.x(), vec.y(), vec.z());
}

template<class T>
Matrix4d rotate3_rpy(const Vector3<T> &vec){
	return rotate3_rpy(vec.x(), vec.y(), vec.z());
}

template<class T>
Matrix4d rotate3_ypr(const Vector3<T> &vec){
	return rotate3_ypr(vec.x(), vec.y(), vec.z());
}

template<class T>
Matrix4d transform3_extrinsic(T x, T y, T z, T rx, T ry, T rz){
	return translate3(x, y, z) * rotate3_xyz(rx, ry, rz);
}

template<class T>
Matrix4d transform3_extrinsic(const Vector3<T> &trans, const Vector3<T> &rot){
	return transform3_extrinsic(trans.x(), trans.y(), trans.z(), rot.x(), rot.y(), rot.z());
}

template<class T>
Matrix4d transform3_intrinsic(T x, T y, T z, T roll, T pitch, T yaw){
	return translate3(x, y, z) * rotate3_rpy(roll, pitch, yaw);
}

template<class T>
Matrix4d transform3_intrinsic(const Vector3<T> &trans, const Vector3<T> &rot){
	return transform3_intrinsic(trans.x(), trans.y(), trans.z(), rot.x(), rot.y(), rot.z());
}

template<class T>
Matrix3d rotation3_axis(const Vector3<T> &axis, const T &angle){
	const auto K = axis.normalized().crossed();
	const auto rot = Matrix3<T>::I() + std::sin(angle) * K + (1 - std::cos(angle)) * (K * K);
	return rot;
}

template<class T>
Matrix3d rotation3_axis(const Vector3<T> &axis_angle){
	const auto phi = axis_angle.norm();
	if(phi < 1e-6){
		return Matrix3<T>::I();
	}
	const auto axis = axis_angle / phi;
	return rotation3_axis<T>(axis, phi);
}


template<class T>
Matrix4d rotate3_axis(const Vector3<T> &axis, const T &angle){
	const auto rot = rotation3_axis(axis, angle);

	Matrix4d result;
	result.set(0, 0, rot);
	result(3, 3) = 1;
	return result;
}

template<class T>
Matrix4d rotate3_axis(const Vector3<T> &axis_angle){
	const auto phi = axis_angle.norm();
	if(phi < 1e-6){
		return Matrix4<T>::I();
	}
	const auto axis = axis_angle / phi;
	return rotate3_axis(axis, phi);
}

template<class T>
Matrix4<T> rotate3_quaternion(T q_0, T q_1, T q_2, T q_3){
	const Matrix4<T> Q = {
		 q_0*q_0+q_1*q_1-q_2*q_2-q_3*q_3, -2*q_0*q_3+2*q_1*q_2,              2*q_0*q_2+2*q_1*q_3,             0,
		 2*q_0*q_3+2*q_1*q_2,              q_0*q_0-q_1*q_1+q_2*q_2-q_3*q_3, -2*q_0*q_1+2*q_2*q_3,             0,
		-2*q_0*q_2+2*q_1*q_3,              2*q_0*q_1+2*q_2*q_3,              q_0*q_0-q_1*q_1-q_2*q_2+q_3*q_3, 0,
		 0,                                0,                                0,                               1,
	};
	return Q;
}

template<class T>
Matrix4<T> rotate3_quaternion(const Vector4<T> &quat){
	const auto quat_ = quat.normalized();
	return rotate3_quaternion(quat_[0], quat_[1], quat_[2], quat_[3]);
}

template<class T>
Matrix4<T> rotate3_quaternion(const Vector3<T> &quat){
	const auto q123_norm = quat.squared_norm();
	if(q123_norm > 1){
		throw std::logic_error("not a union quaternion");
	}
	const auto q_0 = std::sqrt(1 - q123_norm);
	return rotate3_quaternion(q_0, quat[0], quat[1], quat[2]);
}

template<class T>
Vector3<T> get_axis(const Matrix3<T> &rotation){
	Matrix3<T> H = rotation - Matrix3<T>::I();
	if(H.norm() < 1e-6){
		return Vector3<T>();
	}
	H.gauss_jordan_upper();
	const Vector3<T> c_1 = H.template get<1, 3>(0, 0).transposed();
	const Vector3<T> c_2 = H.template get<1, 3>(1, 0).transposed();
	const Vector3<T> axis = c_1.crossed(c_2).normalized();

	const Matrix3<T> N = axis.crossed();
	const auto sin = -(N * rotation).trace() / 2;
	const auto cos = (rotation.trace() - 1) / 2;
	const auto rot = std::atan2(sin, cos);

	return rot * axis;
}

template<class T>
Vector3<T> get_axis(const Matrix4<T> &transform){
	return get_axis(transform.template get<3, 3>(0, 0));
}

template<class T>
Matrix3<T> rotation_interpolated(const Matrix3<T> &first, const Matrix3<T> &second, double t){
	const auto relative = first.inverse() * second;
	const auto axis_angle = get_axis(relative);
	const auto R = rotation3_axis<T>(t * axis_angle);
	return first * R;
}

template<class T>
Matrix4<T> transform_interpolated(const Matrix4<T> &first, const Matrix4<T> &second, double t){
	const auto rot = rotation_interpolated(first.template get<3, 3>(0, 0), second.template get<3, 3>(0, 0), t);
	const auto trans = first.template get<3, 1>(0, 3).linear_interpolated(second.template get<3, 1>(0, 3), t);
	Matrix4<T> result;
	result.set(0, 0, rot);
	result.set(0, 3, trans);
	result(3, 3) = 1;
	return result;
}


}
