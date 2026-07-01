#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <stdexcept>

#include <math/numbers.h>


namespace math{

// Workaround to allow partial template specialization

template<class T, size_t Rows, size_t Cols>
class Matrix;

template<class T, size_t Rows> T determinant(const Matrix<T, Rows, Rows> &A);
template<class T> T determinant(const Matrix<T, 0, 0> &A);
template<class T> T determinant(const Matrix<T, 1, 1> &A);
template<class T> T determinant(const Matrix<T, 2, 2> &A);
template<class T> T determinant(const Matrix<T, 3, 3> &A);
template<class T> T determinant(const Matrix<T, 4, 4> &A);

template<class T, size_t Rows> Matrix<T, Rows, Rows> inverse(const Matrix<T, Rows, Rows> &A);
template<class T> Matrix<T, 0, 0> inverse(const Matrix<T, 0, 0> &A);
template<class T> Matrix<T, 1, 1> inverse(const Matrix<T, 1, 1> &A);
template<class T> Matrix<T, 2, 2> inverse(const Matrix<T, 2, 2> &A);
template<class T> Matrix<T, 3, 3> inverse(const Matrix<T, 3, 3> &A);
template<class T> Matrix<T, 4, 4> inverse(const Matrix<T, 4, 4> &A);



template<class T, size_t Rows, size_t Cols>
class Matrix{
public:
	static const T zero;
	static const T one;

	static Matrix I(){
		static_assert(Rows == Cols, "Rows != Cols");
		Matrix result;
		for(size_t i=0; i<Rows; i++){
			result(i, i) = one;
		}
		return result;
	}

	Matrix(){
		fill(zero);
	}

	Matrix(const Matrix &other) = default;

	template<class U>
	Matrix(const Matrix<U, Rows, Cols> &other){
		*this = other;
	}

	Matrix(Matrix &&other) = default;

	Matrix(const std::initializer_list<T> &list){
		const size_t size = list.size();
		if(size == 0){
			// No initialization
		}else if(size != Rows*Cols){
			throw std::logic_error("wrong initializer list size");
		}else{
			// list in human readable form (row major)
			size_t i=0;
			for(const auto &entry : list){
				(*this)(i/Cols, i%Cols) = entry;
				i++;
			}
		}
	}

	Matrix(const std::array<T, Rows*Cols> &data_):
		data(data_)
	{}

	Matrix(std::array<T, Rows*Cols> &&data_):
		data(std::move(data_))
	{}

	Matrix(const std::vector<T> &data_){
		if(data_.size() != data.size()){
			throw std::logic_error("Wrong data size");
		}
		std::copy(data_.begin(), data_.end(), data.begin());
	}

	Matrix &operator=(const Matrix &other) = default;

	template<class U>
	Matrix &operator=(const Matrix<U, Rows, Cols> &other){
		for(size_t i=0; i<Rows*Cols; i++){
			(*this)[i] = static_cast<T>(other[i]);
		}
		return *this;
	}

	Matrix &operator=(Matrix &&other) = default;

	T &operator()(size_t i, size_t j){
		if(i >= Rows || j >= Cols){
			throw std::logic_error("index out of range");
		}
		return data[i + j*Rows];
	}

	const T &operator()(size_t i, size_t j) const{
		if(i >= Rows || j >= Cols){
			throw std::logic_error("index out of range");
		}
		return data[i + j*Rows];
	}

	T &operator[](size_t i){
		if(i >= Rows*Cols){
			throw std::logic_error("index out of range");
		}
		return data[i];
	}

	const T &operator[](size_t i) const{
		if(i >= Rows*Cols){
			throw std::logic_error("index out of range");
		}
		return data[i];
	}

	template<size_t N, size_t M>
	Matrix<T, N, M> get(size_t i, size_t j) const{
		if(i+N > Rows || j+M > Cols){
			throw std::logic_error("index out of range");
		}
		Matrix<T, N, M> result;
		for(size_t y=0; y<N; y++){
			for(size_t x=0; x<M; x++){
				result(y, x) = (*this)(i+y, j+x);
			}
		}
		return result;
	}

	Matrix<T, 1, Cols> get_row(size_t i) const{
		return get<1, Cols>(i, 0);
	}

	Matrix<T, Rows, 1> get_column(size_t j) const{
		return get<Rows, 1>(0, j);
	}

	Matrix<T, Rows-1, Cols> cancel_row(size_t row) const{
		static_assert(Rows > 0, "Rows <= 0");
		Matrix<T, Rows-1, Cols> result = {};
		size_t i_this = 0;
		for(size_t i=0; i<Rows-1; i++){
			if(i_this == row){
				i_this++;
			}
			for(size_t j=0; j<Cols; j++){
				result(i, j) = (*this)(i_this, j);
			}
			i_this++;
		}
		return result;
	}

	Matrix<T, Rows, Cols-1> cancel_column(size_t column) const{
		static_assert(Cols > 0, "Cols <= 0");
		Matrix<T, Rows, Cols-1> result = {};
		for(size_t i=0; i<Rows; i++){
			size_t j_this = 0;
			for(size_t j=0; j<Cols-1; j++){
				if(j_this == column){
					j_this++;
				}
				result(i, j) = (*this)(i, j_this);
				j_this++;
			}
		}
		return result;
	}

	Matrix<T, Rows-1, Cols-1> cancel(size_t row, size_t column) const{
		static_assert(Rows > 0, "Rows <= 0");
		static_assert(Cols > 0, "Cols <= 0");
		Matrix<T, Rows-1, Cols-1> result = {};
		size_t i_this = 0;
		for(size_t i=0; i<Rows-1; i++){
			if(i_this == row){
				i_this++;
			}
			size_t j_this = 0;
			for(size_t j=0; j<Cols-1; j++){
				if(j_this == column){
					j_this++;
				}
				result(i, j) = (*this)(i_this, j_this);
				j_this++;
			}
			i_this++;
		}
		return result;
	}

	void swap_rows(size_t i_1, size_t i_2){
		const auto tmp = get_row(i_1);
		set(i_1, 0, get_row(i_2));
		set(i_2, 0, tmp);
	}

	void swap_columns(size_t j_1, size_t j_2){
		const auto tmp = get_column(j_1);
		set(0, j_1, get_column(j_2));
		set(0, j_2, tmp);
	}

	template<size_t N, size_t M>
	void set(size_t i, size_t j, const Matrix<T, N, M> &other){
		if(i+N > Rows || j+M > Cols){
			throw std::logic_error("index out of range");
		}
		for(size_t y=0; y<N; y++){
			for(size_t x=0; x<M; x++){
				(*this)(i+y, j+x) = other(y, x);
			}
		}
	}

	void fill(const T &value){
		data.fill(value);
	}

	T squared_norm() const{
		T result = zero;
		for(const auto &entry : data){
			result += entry*entry;
		}
		return result;
	}

	template<class F=typename std::conditional<std::is_floating_point<T>::value, T, double>::type>
	F norm() const{
		return std::sqrt(squared_norm());
	}

	template<class F=typename std::conditional<std::is_floating_point<T>::value, T, double>::type>
	F normalize(){
		const auto n = norm<F>();
		if(n == 0){
			throw std::runtime_error("matrix is zero");
		}
		*this /= n;
		return n;
	}

	template<class F=T>
	Matrix<F, Rows, Cols> normalized(F &norm) const{
		Matrix<F, Rows, Cols> result = *this;
		norm = result.template normalize<F>();
		return result;
	}

	template<class F=T>
	Matrix<F, Rows, Cols> normalized() const{
		F tmp;
		return normalized<F>(tmp);
	}

	T trace() const{
		T result = zero;
		for(size_t i=0; i<std::min(Rows, Cols); i++){
			result += (*this)(i, i);
		}
		return result;
	}

	T determinant() const{
		return math::determinant(*this);
	}

	template<class F=typename std::conditional<std::is_floating_point<T>::value, T, double>::type>
	F inv_determinant() const{
		const auto det = determinant();
		if(std::abs(det - zero) < 1e-6){
			throw std::logic_error("matrix is singular");
		}
		return static_cast<F>(one) / det;
	}

	bool is_invertible() const{
		return (std::abs(determinant() - zero) >= 1e-6);
	}

	Matrix inverse() const{
		return math::inverse(*this);
	}

	void invert(){
		*this = inverse();
	}

	Matrix<T, Rows+1, Cols> extended() const{
		Matrix<T, Rows+1, Cols> result = {};
		result.set(0, 0, *this);
		for(size_t j=0; j<Cols; j++){
			result(Rows, j) = one;
		}
		return result;
	}

	Matrix<T, Rows-1, Cols> projected() const{
		static_assert(Rows > 0, "Rows <= 0");
		Matrix<T, Rows-1, Cols> result = {};
		for(size_t j=0; j<Cols; j++){
			result.set(0, j, get<Rows-1, 1>(0, j) / (*this)(Rows-1, j));
		}
		return result;
	}

	Matrix<T, Cols, Rows> transposed() const{
		Matrix<T, Cols, Rows> result = {};
		for(size_t i=0; i<Rows; i++){
			for(size_t j=0; j<Cols; j++){
				result(j, i) = (*this)(i, j);
			}
		}
		return result;
	}

	void dot(const Matrix &other){
		for(size_t i=0; i<Rows*Cols; i++){
			data[i] *= other.data[i];
		}
	}

	Matrix dotted(const Matrix &other) const{
		Matrix result = *this;
		result.dot(other);
		return result;
	}

	void cross(const Matrix &other){
		static_assert(Rows==3, "Only available in 3 dimensions");
		for(size_t j=0; j<Cols; j++){
			const Matrix<T, Rows, 1> col = {
				(*this)(1, j)*other(2, j) - (*this)(2, j)*other(1, j),
				(*this)(2, j)*other(0, j) - (*this)(0, j)*other(2, j),
				(*this)(0, j)*other(1, j) - (*this)(1, j)*other(0, j),
			};
			set(0, j, col);
		}
	}

	Matrix crossed(const Matrix &other) const{
		Matrix result = *this;
		result.cross(other);
		return result;
	}

	Matrix<T, Rows, Rows> crossed() const{
		static_assert(Rows==3, "Only available in 3 dimensions");
		static_assert(Cols==1, "Only available for vectors");
		const Matrix<T, Rows, Rows> result = {
			0, -(*this)[2], (*this)[1],
			(*this)[2], 0, -(*this)[0],
			-(*this)[1], (*this)[0], 0,
		};
		return result;
	}

	template<class U=T>
	typename std::enable_if<std::is_integral<U>::value, T>::type elimination_target(const T &value, const T &eraser) const{
		return lcm(value, eraser);
	}

	template<class U=T>
	typename std::enable_if<std::is_floating_point<U>::value, T>::type elimination_target(const T &value, const T &/*eraser*/) const{
		return value;
	}

	T eliminate(size_t row, size_t column, size_t with_row){
		return eliminate(row, column, get_row(with_row));
	}

	T eliminate(size_t row, size_t column, const Matrix<T, 1, Cols> &with_row){
		const auto value = (*this)(row, column);
		if(value == zero){
			return one;
		}
		const auto eraser = with_row(0, column);
		if(eraser == zero){
			throw std::logic_error("no elimination possible");
		}
		const auto m = elimination_target(value, eraser);
		const auto row_factor = m / value;
		const auto eraser_factor = m / eraser;
		set(row, 0, row_factor * get_row(row) - eraser_factor * with_row);
		return row_factor;
	}

	// Transforms the matrix into upper triangular form by linear transformations.
	// Returns the factor of the effect on the determinant.
	T gauss_jordan_upper(){
		T factor = one;
		size_t step_row = 0;
		size_t step_col = 0;
		while(step_row < Rows && step_col < Cols){
			size_t pivot_row = step_row;
			T pivot = (*this)(pivot_row, step_col);
			for(size_t i=pivot_row+1; i<Rows; i++){
				const T candidate = (*this)(i, step_col);
				if(std::abs(candidate) > std::abs(pivot)){
					pivot_row = i;
					pivot = candidate;
				}
			}
			if(pivot == zero){
				step_col++;
				continue;
			}
			if(pivot_row > step_row){
				swap_rows(pivot_row, step_row);
				factor = -factor;
			}
			const auto row = get_row(step_row);
			for(size_t lower=step_row+1; lower<Rows; lower++){
				factor *= eliminate(lower, step_col, row);
			}
			step_row++;
			step_col++;
		}
		return factor;
	}

	// Completes gaussian eliminination and transforms the matrix into row-echelon form.
	// Assumes that the matrix is already in upper triangular form, like after calling gauss_jordan_upper().
	// Returns the factor of the effect on the determinant.
	T gauss_jordan_lower(){
		T factor = one;
		size_t step_row = 0;
		size_t step_col = 0;
		while(step_row < Rows && step_col < Cols){
			if((*this)(step_row, step_col) == zero){
				step_col++;
				continue;
			}
			const auto row = get_row(step_row);
			for(size_t upper=0; upper<step_row; upper++){
				factor *= eliminate(upper, step_col, row);
			}
			step_row++;
			step_col++;
		}
		return factor;
	}

	void gauss_jordan(){
		gauss_jordan_upper();
		gauss_jordan_lower();
	}

	Matrix &operator+=(const Matrix &other){
		for(size_t i=0; i<Rows*Cols; i++){
			data[i] += other.data[i];
		}
		return *this;
	}

	Matrix operator+(const Matrix &other) const{
		Matrix result = *this;
		result += other;
		return result;
	}

	Matrix &operator-=(const Matrix &other){
		for(size_t i=0; i<Rows*Cols; i++){
			data[i] -= other.data[i];
		}
		return *this;
	}

	Matrix operator-(const Matrix &other) const{
		Matrix result = *this;
		result -= other;
		return result;
	}

	Matrix &operator*=(const T &scalar){
		for(auto &entry : data){
			entry *= scalar;
		}
		return *this;
	}

	Matrix operator*(const T &scalar) const{
		Matrix result = *this;
		result *= scalar;
		return result;
	}

	template<size_t M>
	Matrix<T, Rows, M> operator*(const Matrix<T, Cols, M> &other) const{
		Matrix<T, Rows, M> result;
		for(size_t i=0; i<Rows; i++){
			for(size_t j=0; j<M; j++){
				auto &cell = result(i, j);
				for(size_t k=0; k<Cols; k++){
					cell += (*this)(i, k) * other(k, j);
				}
			}
		}
		return result;
	}

	Matrix &operator/=(const T &scalar){
		for(auto &entry : data){
			entry /= scalar;
		}
		return *this;
	}

	Matrix operator/(const T &scalar) const{
		Matrix result = *this;
		result /= scalar;
		return result;
	}

	bool operator==(const Matrix &other) const{
		return data == other.data;
	}

	bool operator!=(const Matrix &other) const{
		return data != other.data;
	}

private:
	// data stored in column major
	std::array<T, Rows*Cols> data;
};

template<class T, size_t Rows, size_t Cols>
const T Matrix<T, Rows, Cols>::zero = static_cast<T>(0);

template<class T, size_t Rows, size_t Cols>
const T Matrix<T, Rows, Cols>::one = static_cast<T>(1);


template<class T>
T determinant(const Matrix<T, 0, 0> &A){
	return A.one;
}

template<class T>
T determinant(const Matrix<T, 1, 1> &A){
	return A(0, 0);
}

template<class T>
T determinant(const Matrix<T, 2, 2> &A){
	return A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
}

template<class T>
T determinant(const Matrix<T, 3, 3> &A){
	return A.zero
		+ A(0, 0) * A(1, 1) * A(2, 2)
		+ A(0, 1) * A(1, 2) * A(2, 0)
		+ A(0, 2) * A(1, 0) * A(2, 1)
		- A(0, 0) * A(1, 2) * A(2, 1)
		- A(0, 1) * A(1, 0) * A(2, 2)
		- A(0, 2) * A(1, 1) * A(2, 0)
	;
}

template<class T>
T determinant(const Matrix<T, 4, 4> &A){
	return A.zero
		+ A(0, 0) * (0
			+ A(1, 1) * A(2, 2) * A(3, 3)
			+ A(1, 2) * A(2, 3) * A(3, 1)
			+ A(1, 3) * A(2, 1) * A(3, 2)
			- A(1, 1) * A(2, 3) * A(3, 2)
			- A(1, 2) * A(2, 1) * A(3, 3)
			- A(1, 3) * A(2, 2) * A(3, 1)
		)
		- A(0, 1) * (0
			+ A(1, 0) * A(2, 2) * A(3, 3)
			+ A(1, 2) * A(2, 3) * A(3, 0)
			+ A(1, 3) * A(2, 0) * A(3, 2)
			- A(1, 0) * A(2, 3) * A(3, 2)
			- A(1, 2) * A(2, 0) * A(3, 3)
			- A(1, 3) * A(2, 2) * A(3, 0)
		)
		+ A(0, 2) * (0
			+ A(1, 0) * A(2, 1) * A(3, 3)
			+ A(1, 1) * A(2, 3) * A(3, 0)
			+ A(1, 3) * A(2, 0) * A(3, 1)
			- A(1, 0) * A(2, 3) * A(3, 1)
			- A(1, 1) * A(2, 0) * A(3, 3)
			- A(1, 3) * A(2, 1) * A(3, 0)

		)
		- A(0, 3) * (0
			+ A(1, 0) * A(2, 1) * A(3, 2)
			+ A(1, 1) * A(2, 2) * A(3, 0)
			+ A(1, 2) * A(2, 0) * A(3, 1)
			- A(1, 0) * A(2, 2) * A(3, 1)
			- A(1, 1) * A(2, 0) * A(3, 2)
			- A(1, 2) * A(2, 1) * A(3, 0)
		)
	;

}

template<class T, size_t Rows>
T determinant(const Matrix<T, Rows, Rows> &A){
	Matrix<T, Rows, Rows> tmp = A;
	const T factor = tmp.gauss_jordan_upper();
	T det = A.one;
	for(size_t i=0; i<Rows; i++){
		det *= tmp(i, i);
	}
	return det / factor;

	/*
	// Laplace expansion along the first column (very inefficient)
	T result = A.zero;
	const auto canceled = A.cancel_column(0);
	for(size_t i=0; i<Rows; i++){
		const auto summand = A(i, 0) * canceled.cancel_row(i).determinant();
		if(i % 2 == 0){
			result += summand;
		}else{
			result -= summand;
		}
	}
	return result;
	*/
}

template<class T>
Matrix<T, 0, 0> inverse(const Matrix<T, 0, 0> &A){
	return A;
}

template<class T>
Matrix<T, 1, 1> inverse(const Matrix<T, 1, 1> &A){
	const auto inv_determinant = A.inv_determinant();
	return {
		inv_determinant,
	};
}

template<class T>
Matrix<T, 2, 2> inverse(const Matrix<T, 2, 2> &A){
	const auto inv_determinant = A.inv_determinant();
	const Matrix<T, 2, 2> adj = {
		 A(1, 1), -A(0, 1),
		-A(1, 0),  A(0, 0),
	};
	return inv_determinant * adj;
}

template<class T>
Matrix<T, 3, 3> inverse(const Matrix<T, 3, 3> &A){
	const auto inv_determinant = A.inv_determinant();
	const Matrix<T, 3, 3> adj = {
		A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1),  A(0, 2)*A(2, 1) - A(0, 1)*A(2, 2),  A(0, 1)*A(1, 2) - A(0, 2)*A(1, 1),
		A(1, 2)*A(2, 0) - A(1, 0)*A(2, 2),  A(0, 0)*A(2, 2) - A(0, 2)*A(2, 0),  A(0, 2)*A(1, 0) - A(0, 0)*A(1, 2),
		A(1, 0)*A(2, 1) - A(1, 1)*A(2, 0),  A(0, 1)*A(2, 0) - A(0, 0)*A(2, 1),  A(0, 0)*A(1, 1) - A(0, 1)*A(1, 0),
	};
	return inv_determinant * adj;
}

template<class T>
Matrix<T, 4, 4> inverse(const Matrix<T, 4, 4> &A){
	const auto inv_determinant = A.inv_determinant();
	const Matrix<T, 4, 4> adj = {
		 A(1, 1)*A(2, 2)*A(3, 3) + A(1, 2)*A(2, 3)*A(3, 1) + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1) - A(1, 2)*A(2, 1)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2),
		-A(0, 1)*A(2, 2)*A(3, 3) - A(0, 2)*A(2, 3)*A(3, 1) - A(0, 3)*A(2, 1)*A(3, 2) + A(0, 3)*A(2, 2)*A(3, 1) + A(0, 2)*A(2, 1)*A(3, 3) + A(0, 1)*A(2, 3)*A(3, 2),
		 A(0, 1)*A(1, 2)*A(3, 3) + A(0, 2)*A(1, 3)*A(3, 1) + A(0, 3)*A(1, 1)*A(3, 2) - A(0, 3)*A(1, 2)*A(3, 1) - A(0, 2)*A(1, 1)*A(3, 3) - A(0, 1)*A(1, 3)*A(3, 2),
		-A(0, 1)*A(1, 2)*A(2, 3) - A(0, 2)*A(1, 3)*A(2, 1) - A(0, 3)*A(1, 1)*A(2, 2) + A(0, 3)*A(1, 2)*A(2, 1) + A(0, 2)*A(1, 1)*A(2, 3) + A(0, 1)*A(1, 3)*A(2, 2),

		-A(1, 0)*A(2, 2)*A(3, 3) - A(1, 2)*A(2, 3)*A(3, 0) - A(1, 3)*A(2, 0)*A(3, 2) + A(1, 3)*A(2, 2)*A(3, 0) + A(1, 2)*A(2, 0)*A(3, 3) + A(1, 0)*A(2, 3)*A(3, 2),
		 A(0, 0)*A(2, 2)*A(3, 3) + A(0, 2)*A(2, 3)*A(3, 0) + A(0, 3)*A(2, 0)*A(3, 2) - A(0, 3)*A(2, 2)*A(3, 0) - A(0, 2)*A(2, 0)*A(3, 3) - A(0, 0)*A(2, 3)*A(3, 2),
		-A(0, 0)*A(1, 2)*A(3, 3) - A(0, 2)*A(1, 3)*A(3, 0) - A(0, 3)*A(1, 0)*A(3, 2) + A(0, 3)*A(1, 2)*A(3, 0) + A(0, 2)*A(1, 0)*A(3, 3) + A(0, 0)*A(1, 3)*A(3, 2),
		 A(0, 0)*A(1, 2)*A(2, 3) + A(0, 2)*A(1, 3)*A(2, 0) + A(0, 3)*A(1, 0)*A(2, 2) - A(0, 3)*A(1, 2)*A(2, 0) - A(0, 2)*A(1, 0)*A(2, 3) - A(0, 0)*A(1, 3)*A(2, 2),

		 A(1, 0)*A(2, 1)*A(3, 3) + A(1, 1)*A(2, 3)*A(3, 0) + A(1, 3)*A(2, 0)*A(3, 1) - A(1, 3)*A(2, 1)*A(3, 0) - A(1, 1)*A(2, 0)*A(3, 3) - A(1, 0)*A(2, 3)*A(3, 1),
		-A(0, 0)*A(2, 1)*A(3, 3) - A(0, 1)*A(2, 3)*A(3, 0) - A(0, 3)*A(2, 0)*A(3, 1) + A(0, 3)*A(2, 1)*A(3, 0) + A(0, 1)*A(2, 0)*A(3, 3) + A(0, 0)*A(2, 3)*A(3, 1),
		 A(0, 0)*A(1, 1)*A(3, 3) + A(0, 1)*A(1, 3)*A(3, 0) + A(0, 3)*A(1, 0)*A(3, 1) - A(0, 3)*A(1, 1)*A(3, 0) - A(0, 1)*A(1, 0)*A(3, 3) - A(0, 0)*A(1, 3)*A(3, 1),
		-A(0, 0)*A(1, 1)*A(2, 3) - A(0, 1)*A(1, 3)*A(2, 0) - A(0, 3)*A(1, 0)*A(2, 1) + A(0, 3)*A(1, 1)*A(2, 0) + A(0, 1)*A(1, 0)*A(2, 3) + A(0, 0)*A(1, 3)*A(2, 1),

		-A(1, 0)*A(2, 1)*A(3, 2) - A(1, 1)*A(2, 2)*A(3, 0) - A(1, 2)*A(2, 0)*A(3, 1) + A(1, 2)*A(2, 1)*A(3, 0) + A(1, 1)*A(2, 0)*A(3, 2) + A(1, 0)*A(2, 2)*A(3, 1),
		 A(0, 0)*A(2, 1)*A(3, 2) + A(0, 1)*A(2, 2)*A(3, 0) + A(0, 2)*A(2, 0)*A(3, 1) - A(0, 2)*A(2, 1)*A(3, 0) - A(0, 1)*A(2, 0)*A(3, 2) - A(0, 0)*A(2, 2)*A(3, 1),
		-A(0, 0)*A(1, 1)*A(3, 2) - A(0, 1)*A(1, 2)*A(3, 0) - A(0, 2)*A(1, 0)*A(3, 1) + A(0, 2)*A(1, 1)*A(3, 0) + A(0, 1)*A(1, 0)*A(3, 2) + A(0, 0)*A(1, 2)*A(3, 1),
		 A(0, 0)*A(1, 1)*A(2, 2) + A(0, 1)*A(1, 2)*A(2, 0) + A(0, 2)*A(1, 0)*A(2, 1) - A(0, 2)*A(1, 1)*A(2, 0) - A(0, 1)*A(1, 0)*A(2, 2) - A(0, 0)*A(1, 2)*A(2, 1),
	};
	return inv_determinant * adj;
}

template<class T, size_t Rows>
Matrix<T, Rows, Rows> inverse(const Matrix<T, Rows, Rows> &A){
	const auto id = Matrix<T, Rows, Rows>::I();
	Matrix<T, Rows, 2*Rows> H = {};
	H.set(0, 0, A);
	H.set(0, Rows, id);
	H.gauss_jordan_upper();
	for(size_t i=0; i<Rows; i++){
		if(H(i, i) == A.zero){
			throw std::logic_error("matrix is singular");
		}
	}
	H.gauss_jordan_lower();
	for(size_t i=0; i<Rows; i++){
		H.set(i, 0, H.get_row(i) / H(i, i));
	}
	return H.template get<Rows, Rows>(0, Rows);
}


template<class T, size_t Rows, size_t Cols>
Matrix<T, Rows, Cols> operator*(const T &scalar, const Matrix<T, Rows, Cols> &matrix){
	return matrix * scalar;
}


template<class T> using Matrix2 = Matrix<T, 2, 2>;
template<class T> using Matrix3 = Matrix<T, 3, 3>;
template<class T> using Matrix4 = Matrix<T, 4, 4>;

using Matrix2i = Matrix2<int>;
using Matrix2f = Matrix2<float>;
using Matrix2d = Matrix2<double>;
using Matrix3i = Matrix3<int>;
using Matrix3f = Matrix3<float>;
using Matrix3d = Matrix3<double>;
using Matrix4i = Matrix4<int>;
using Matrix4f = Matrix4<float>;
using Matrix4d = Matrix4<double>;

}
