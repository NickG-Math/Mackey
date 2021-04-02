#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen won't work otherwise.
#endif

#ifdef _OPENMP
#define MACKEY_USE_OPEN_MP
#endif

//Eigen directives
//#define EIGEN_USE_MKL_ALL
//#define EIGEN_NO_DEBUG //disable Eigen assertions
//#define NDEBUG
#include "Coefficients/Z_n.hpp"
#include "Groups/C4.hpp"
#include <iostream>
#include <array>
#include <chrono>
#include "C4Verify.ipp"

typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::Eigen<short, -1, -1> dense_diff_t;
typedef Eigen::SparseMatrix<short> sparse_diff_t;
typedef C4<rank_t, dense_diff_t> dense_group_t;
typedef C4<rank_t, sparse_diff_t> sparse_group_t;

int main()
{
	constexpr bool user_facing = 1;

	constexpr bool verify = 1;
	constexpr bool verify_deep = 1;

	constexpr bool performance_test = 1;

	if constexpr (user_facing)
	{
		basic_C4_add_mul_Z_coeff_test<dense_group_t>(1, 0, 0, 0);
		return 0;
	}
	
	if constexpr (verify)
	{
		basic_C4_add_mul_Z_coeff_test<dense_group_t>(0, 1, !verify_deep, 0);
	}

	if constexpr (performance_test)
	{
			std::cout << "\n PERFORMANCE DENSE" << std::endl;
			basic_C4_add_mul_Z_coeff_test<dense_group_t>(0, 0, 0, 1);
			std::cout << "\n PERFORMANCE SPARSE" << std::endl;
			basic_C4_add_mul_Z_coeff_test<sparse_group_t>(0, 0, 0, 1);
	}

	return 0;
}
