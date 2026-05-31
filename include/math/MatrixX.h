#pragma once

#include <vector>
#include <stdexcept>


namespace math{


template<class T>
class MatrixX{
public:
	MatrixX(size_t rows, size_t cols):
		Rows(rows),
		Cols(cols)
	{
		data.resize(rows*cols);
	}

	size_t rows() const{
		return Rows;
	}

	size_t cols() const{
		return Cols;
	}

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

	std::vector<T> get_data() const{
		return data;
	}

	MatrixX<T> transposed() const{
		MatrixX<T> result(Cols, Rows);
		for(size_t i=0; i<Rows; i++){
			for(size_t j=0; j<Cols; j++){
				result(j, i) = (*this)(i, j);
			}
		}
		return result;
	}

	MatrixX &operator=(const std::initializer_list<T> &list){
		const size_t size = list.size();
		if(size == 0){
			// no initialization
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
		return *this;
	}

	MatrixX &operator+=(const MatrixX &other){
		if(other.rows() != Rows || other.cols() != Cols){
			throw std::logic_error("Dimensions do not match");
		}
		for(size_t i=0; i<Rows*Cols; i++){
			data[i] += other.data[i];
		}
		return *this;
	}

	MatrixX operator+(const MatrixX &other) const{
		MatrixX result = *this;
		result += other;
		return result;
	}

	MatrixX &operator-=(const MatrixX &other){
		if(other.rows() != Rows || other.cols() != Cols){
			throw std::logic_error("Dimensions do not match");
		}
		for(size_t i=0; i<Rows*Cols; i++){
			data[i] -= other.data[i];
		}
		return *this;
	}

	MatrixX operator-(const MatrixX &other) const{
		MatrixX result = *this;
		result -= other;
		return result;
	}

	MatrixX &operator*=(const T &scalar){
		for(auto &entry : data){
			entry *= scalar;
		}
		return *this;
	}

	MatrixX operator*(const T &scalar) const{
		MatrixX result = *this;
		result *= scalar;
		return result;
	}

	MatrixX<T> operator*(const MatrixX<T> &other) const{
		const size_t M = other.cols();
		if(M != Rows){
			throw std::logic_error("Dimensions do not match");
		}
		MatrixX<T> result(Rows, M);
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

	MatrixX &operator/=(const T &scalar){
		for(auto &entry : data){
			entry /= scalar;
		}
		return *this;
	}

	MatrixX operator/(const T &scalar) const{
		MatrixX result = *this;
		result /= scalar;
		return result;
	}

	bool operator==(const MatrixX &other) const{
		if(other.rows() != Rows || other.cols() != Cols){
			return false;
		}
		return data == other.data;
	}

	bool operator!=(const MatrixX &other) const{
		if(other.rows() != Rows || other.cols() != Cols){
			return false;
		}
		return data != other.data;
	}

private:
	const size_t Rows;
	const size_t Cols;
	// data stored in column major
	std::vector<T> data;
};



}
