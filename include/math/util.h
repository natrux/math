#pragma once

#include <math/Matrix.h>

#include <iostream>


namespace math{


template<class T, size_t N, size_t M>
void print_matrix(const Matrix<T, N, M> &matrix){
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


}

