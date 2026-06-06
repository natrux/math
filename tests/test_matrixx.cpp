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


