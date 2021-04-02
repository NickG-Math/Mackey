#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen won't work otherwise.
#endif
#define MACKEY_USE_OPEN_MP
#include "Groups/C4.hpp"
#include "impl/Test_C4.ipp"
using namespace mackey;
int main()
{
	test::basic_C4_add_mul_Z_coeff_test<C4<Eigen::Matrix<short, 1, -1>, Eigen::SparseMatrix<short>>>(1, 0, 0, 0);
	return 0;
}
