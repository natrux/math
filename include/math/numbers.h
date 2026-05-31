#pragma once


namespace math{


template<class T>
typename std::enable_if<std::is_integral<T>::value, T>::type gcd(T a, T b){
	while(b != 0){
		const auto h = a % b;
		a = b;
		b = h;
	}
	return std::abs(a);
}


template<class T, class... U>
typename std::enable_if<std::is_integral<T>::value, T>::type gcd(T first, T second, T third, U... tail){
	return gcd(gcd(first, second), third, tail...);
}


template<class T>
typename std::enable_if<std::is_integral<T>::value, T>::type lcm(T a, T b){
	const auto g = gcd(a, b);
	const T tmp = (std::abs(a) > std::abs(b)) ? (a / g * b) : (b / g * a);
	return std::abs(tmp);
}


template<class T, class... U>
typename std::enable_if<std::is_integral<T>::value, T>::type lcm(T first, T second, T third, U... tail){
	return lcm(lcm(first, second), third, tail...);
}


}
