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
#include "impl/Test_Space.ipp"

int main()
{
	mackey::test::test_BC4S2<C4<Eigen::Matrix<short, 1, -1>, Eigen::SparseMatrix<Z2>>>();
	return 0;
}