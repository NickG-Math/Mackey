//#define EIGEN_USE_MKL_ALL

#ifdef _MSC_VER
#define _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING //Eigen won't work otherwise.
#endif

#include <iostream>
#include "C2n_Implementation.h"
#include "Factorization.h"

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::Matrix<Z<2>, -1, -1> diff_t;

int main() {
	

	std::vector<std::vector<int>> true_sources = { {0,0,0}, {-3,0,-2} , {-2,0,-1} }; 
	std::vector<std::string> source_names = { "1" , "s", "theta" };

	//Go from S^{-5*sigma-5*lambda} to S^{5*sigma+5*lambda}. The basic irreducibles are chosen together with their names
	Factorization<rank_t, diff_t> F(2,{ -5,-5 }, {5,5}, { {0,1,0},{1,1,0},{0,0,1},{2,0,1} }, { "asigma", "usigma", "alambda", "ulambda" });
	//Factorization<rank_t, diff_t> F(1, { -5,-5 }, { 5,5 }, {{1,1,0},{0,0,1},{2,0,1} }, {"usigma", "alambda", "ulambda" }); //middle level

	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	std::cout << F.generator_names();
	std::ofstream file;
	file.open("graph.dot");
	file << (F.getgraph() << F.generator_names());
	file.close();
	return 0;
}
