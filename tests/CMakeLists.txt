cmake_minimum_required(VERSION 3.22)


find_package(Catch2 3 REQUIRED)


add_executable(test_numbers test_numbers.cpp)
target_link_libraries(test_numbers
	math
	Catch2::Catch2WithMain
)


add_executable(test_matrix test_matrix.cpp)
target_link_libraries(test_matrix
	math
	Catch2::Catch2WithMain
)


add_executable(test_vector test_vector.cpp)
target_link_libraries(test_vector
	math
	Catch2::Catch2WithMain
)


add_executable(test_matrixx test_matrixx.cpp)
target_link_libraries(test_matrixx
	math
	Catch2::Catch2WithMain
)


add_executable(test_geometry test_geometry.cpp)
target_link_libraries(test_geometry
	math
	Catch2::Catch2WithMain
)


add_test(NAME test_numbers COMMAND test_numbers)
add_test(NAME test_matrix COMMAND test_matrix)
add_test(NAME test_vector COMMAND test_vector)
add_test(NAME test_matrixx COMMAND test_matrixx)
add_test(NAME test_geometry COMMAND test_geometry)

