#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen won't work otherwise.
#endif
#define MACKEY_USE_OPEN_MP
#include "Groups/C4.hpp"
#include "impl/Test_Massey.ipp"


int main()
{
	mackey::test::test_Massey<C4<Eigen::Matrix<short, 1, -1>, Eigen::SparseMatrix<short>>>();
	return 0;
}