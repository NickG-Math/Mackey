//Eigen directives
#define EIGEN_USE_MKL_ALL
//#define EIGEN_NO_DEBUG //disable Eigen assertions
//#define NDEBUG

#ifdef _MSC_VER
#define _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING //Eigen won't work otherwise.
#endif
//note: These preprocessor directives must be before any includes.

#include <iostream>
#include <chrono>
#include "Implementation.h"
#include "C4Verify.h" //to verify

typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::Matrix<char, -1, -1> diff_t;

int main() {
	auto begin = std::chrono::high_resolution_clock::now();
	C4Test::C4Multtest<rank_t, diff_t>(13, 13, 13, 13, 1);
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;
}
