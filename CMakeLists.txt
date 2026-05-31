cmake_minimum_required(VERSION 3.22)
project(math CXX)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(
	-pedantic -Wall -Wextra
	-Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization
	-Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations
	-Wmissing-include-dirs -Wold-style-cast -Wredundant-decls -Wshadow
	-Wsign-conversion -Wsign-promo -Wstrict-null-sentinel
	-Wstrict-overflow -Wsuggest-override -Wswitch-default -Wundef

	#-Werror
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


add_executable(test
	test.cpp
)

target_include_directories(test PUBLIC
	include
)

