#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen won't work otherwise.
#endif

#define MACKEY_USE_OPEN_MP
#include "Mackey.hpp"
#include "impl/Test_C4.ipp"
#include "impl/Test_Factorization.ipp"
#include "impl/Test_Massey.ipp"
#include "impl/Test_Space.ipp"

using namespace mackey;

typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::Matrix<short, -1, -1> dZdiff_t;
typedef Eigen::Matrix<Z2, -1, -1> dZ2diff_t;
typedef Eigen::SparseMatrix<Z2> sZ2diff_t;

template <bool Z_coeff, bool dense>
using group_t = C4<rank_t, std::conditional_t<dense, Eigen::Matrix<std::conditional_t<Z_coeff, char, Z2>, -1, -1>, Eigen::SparseMatrix<std::conditional_t<Z_coeff, char, Z2>>>>;

int main()
{
	constexpr bool user_facing = 0;

	constexpr bool verify = 1;
	constexpr bool verify_deep = 1;

	constexpr bool performance_test = 1;

	constexpr bool factorization = 1;
	constexpr bool massey = 1;
	constexpr bool space = 1;

	constexpr bool zcoeff = 1;
	constexpr bool z2coeff = 1;

	constexpr bool dense = 1;
	constexpr bool sparse = 1;

	if constexpr (user_facing)
	{
		mackey::test::basic_C4_add_mul_Z_coeff_test<group_t<1, 1>>(1, 0, 0, 0);
		return 0;
	}
	if constexpr (verify)
	{
		//if constexpr (dense) {
		std::cout << "\n VERIFY DENSE\n";
		mackey::test::basic_C4_add_mul_Z_coeff_test<group_t<1, 1>>(0, 1, !verify_deep, 0);
		//}
		//if constexpr (sparse) {
		std::cout << "\n VERIFY SPARSE\n";
		mackey::test::basic_C4_add_mul_Z_coeff_test<group_t<1, 0>>(0, 1, !verify_deep, 0);
		//}
	}

	if constexpr (performance_test)
	{
		if constexpr (dense)
		{
			std::cout << "\n PERFORMANCE DENSE" << std::endl;
			mackey::test::basic_C4_add_mul_Z_coeff_test<group_t<1, 1>>(0, 0, 0, 1);
		}
		if constexpr (performance_test)
		{
			std::cout << "\n PERFORMANCE SPARSE" << std::endl;
			mackey::test::basic_C4_add_mul_Z_coeff_test<group_t<1, 0>>(0, 0, 0, 1);
		}
	}

	if constexpr (factorization)
	{
		if constexpr (zcoeff)
		{
			if constexpr (dense)
			{
				std::cout << "\n FACTORIZATION Z DENSE\n";
				mackey::test::test_factorization<group_t<1, 1>>(8, 0,0);
			}
			if constexpr (sparse)
			{
				std::cout << "\n FACTORIZATION Z SPARSE\n";
				mackey::test::test_factorization<group_t<1, 0>>(8, 0,0);
			}
		}
		if constexpr (z2coeff)
		{
			if constexpr (dense)
			{
				std::cout << "\n FACTORIZATION Z2 DENSE\n";
				mackey::test::test_factorization<group_t<0, 1>>(4, 0,0);
			}
			if constexpr (sparse)
			{
				std::cout << "\n FACTORIZATION Z2 SPARSE\n";
				mackey::test::test_factorization<group_t<0, 0>>(4, 0,1);
			}
		}
	}

	if constexpr (massey)
	{
		if constexpr (zcoeff)
		{
			if constexpr (dense)
			{
				std::cout << "\n MASSEY Z DENSE\n";
				mackey::test::test_Massey<group_t<1, 1>>();
			}
			if constexpr (sparse)
			{
				std::cout << "\n MASSEY Z SPARSE\n";
				mackey::test::test_Massey<group_t<1, 0>>();
			}
		}
		if constexpr (z2coeff)
		{
			if constexpr (dense)
			{
				std::cout << "\n MASSEY Z2 DENSE\n";
				mackey::test::test_Massey<group_t<0, 1>>();
			}
			if constexpr (sparse)
			{
				std::cout << "\n MASSEY Z2 SPARSE\n";
				mackey::test::test_Massey<group_t<0, 0>>();
			}
		}
	}

	if constexpr (space)
	{
		if constexpr (z2coeff)
		{
			if constexpr (dense)
			{
				std::cout << "\n SPACE DENSE\n";
				mackey::test::test_BC4S2<group_t<0, 1>>();
			}
			if constexpr (sparse)
			{
				std::cout << "\n SPACE SPARSE\n";
				mackey::test::test_BC4S2<group_t<0, 0>>();
			}
		}
	}

	return 0;
}
