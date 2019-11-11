#define EIGEN_USE_MKL_ALL ///<Use the Intel MKL with Eigen
//#define EIGEN_NO_DEBUG ///<Disable Eigen asserations
//#define NDEBUG

//#define ADDITIVE_STRUCTURE 		///<Compute the additive structure
#define MULTIPLICATIVE_STRUCTURE	///<Compute the multiplicative structure

#ifdef ADDITIVE_STRUCTURE
#define MACKEY_NAMES ///<If defined, a list of Mackey functors and their names in Group_Specific_Implementations_Optional is going to be used; make sure it's there.
#endif

#ifdef _MSC_VER
#define _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING ///<Eigen won't work otherwise with MSVC.
#endif
//note: These preprocessor directives must be before any includes.

#include <iostream>
#include <chrono>
#include "Implementation.h"
#include "Optional_Implementation.h"
#include  "C4Verify.h" //to verify

typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::Matrix<char, -1, -1> diff_t;


int main() {
	auto begin = std::chrono::high_resolution_clock::now();
#ifdef ADDITIVE_STRUCTURE
	C4Test::C4MackeyTest<rank_t, diff_t>(-40, 0, 0, 30, 0);
#endif
#ifdef MULTIPLICATIVE_STRUCTURE
	C4Test::C4Multtest<rank_t, diff_t>(13, 13, 13, 13, 0);
#endif
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;
	return 0;
}
