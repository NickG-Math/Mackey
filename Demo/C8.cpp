//Eigen directives
#define EIGEN_USE_MKL_ALL
//#define EIGEN_NO_DEBUG //disable Eigen assertions
//#define NDEBUG

#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen/cereal won't work otherwise.
#endif
//note: These preprocessor directives must be before any includes.

#include <iostream>
#include <fstream>
#include "C2n_Implementation.h"
#include <Mackey/Cerealizer.h>
#include "omp.h"
#include <chrono>


typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::Matrix<short, -1, -1> diff_t;

using namespace Mackey;


int main() {

	std::vector<int> minimum = { -6,-6,-6 };
	std::vector<int> maximum = { 6,6,6 };

	//std::vector<int> minimum(3);
	//std::vector<int> maximum(3);
	//std::cout << "min:";
	//for (int i = 0; i < 3; i++) {
	//	std::cin >> minimum[i];
	//}
	//std::cout << "max:";
	//for (int i = 0; i < 3; i++) {
	//	std::cin >> maximum[i];
	//}
	auto begin = std::chrono::high_resolution_clock::now();

	AdditiveStructure<rank_t, diff_t> A(minimum, maximum);

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;

	//saver(A, "add.bin");

	//AdditiveStructure<rank_t, diff_t> A;
	//loader(A, "add.bin");
	A.identify();

	std::ofstream out("C8.txt", std::ofstream::out);
	A.print_answer(out);
	std::ofstream outer("C8names.txt", std::ofstream::out);
	A.print_unique(outer);
	std::stringstream outer2;
	A.print_unknown(outer2);
	std::cout << outer2.str();
}
