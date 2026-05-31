#pragma once

#include <math/Matrix.h>

namespace math{


template<class T, size_t N> using Vector = Matrix<T, N, 1>;


template<class T>
class Vector2 : public Vector<T, 2>{
public:
	Vector2() : Vector<T, 2>() {}
	Vector2(const Vector<T, 2> &other) : Vector<T, 2>(other) {}
	Vector2(const std::initializer_list<T> &list) : Vector<T, 2>(list) {}

	T &x(){ return (*this)[0]; }
	const T &x() const{ return (*this)[0]; }
	T &y(){ return (*this)[1]; }
	const T &y() const{ return (*this)[1]; }

	template<class F=typename std::conditional<std::is_floating_point<T>::value, T, double>::type>
	F atan2() const{
		return std::atan2(y(), x());
	}
};


template<class T>
class Vector3 : public Vector<T, 3>{
public:
	Vector3() : Vector<T, 3>() {}
	Vector3(const Vector<T, 3> &other) : Vector<T, 3>(other) {}
	Vector3(const std::initializer_list<T> &list) : Vector<T, 3>(list) {}

	T &x(){ return (*this)[0]; }
	const T &x() const{ return (*this)[0]; }
	T &y(){ return (*this)[1]; }
	const T &y() const{ return (*this)[1]; }
	T &z(){ return (*this)[2]; }
	const T &z() const{ return (*this)[2]; }
};


template<class T>
class Vector4 : public Vector<T, 4>{
public:
	Vector4() : Vector<T, 4>() {}
	Vector4(const Vector<T, 4> &other) : Vector<T, 4>(other) {}
	Vector4(const std::initializer_list<T> &list) : Vector<T, 4>(list) {}

	T &x(){ return (*this)[0]; }
	const T &x() const{ return (*this)[0]; }
	T &y(){ return (*this)[1]; }
	const T &y() const{ return (*this)[1]; }
	T &z(){ return (*this)[2]; }
	const T &z() const{ return (*this)[2]; }
	T &w(){ return (*this)[3]; }
	const T &w() const{ return (*this)[3]; }
};


using Vector2i = Vector2<int>;
using Vector2f = Vector2<float>;
using Vector2d = Vector2<double>;
using Vector3i = Vector3<int>;
using Vector3f = Vector3<float>;
using Vector3d = Vector3<double>;
using Vector4i = Vector4<int>;
using Vector4f = Vector4<float>;
using Vector4d = Vector4<double>;

}
