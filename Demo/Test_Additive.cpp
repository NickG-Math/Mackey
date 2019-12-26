//Eigen directives
//#define EIGEN_USE_MKL_ALL
//#define EIGEN_NO_DEBUG //disable Eigen assertions
//#define NDEBUG

#define MACKEY_NAMES

#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen/cereal won't work otherwise.
#endif
//note: These preprocessor directives must be before any includes.

#include <iostream>
#include <chrono>
#include <Mackey/Cerealizer.h>
#include "C4_Implementation.h"
#include "C4Verify.h" //to verify
#include <Mackey/Z_n.h>


typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::Matrix<char, -1, -1> diff_t;


int main() {
	auto begin = std::chrono::high_resolution_clock::now();
	C4Test::C4MackeyTest<rank_t, diff_t>(-30, 30, -30, 30, 0);
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;
}
