//#define EIGEN_USE_MKL_ALL

#include <iostream>
#include "Implementation.h"
#include "Mackey/Factorization.h"

using namespace Mackey;

typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::Matrix<char, -1, -1> diff_t;

int main() {

	std::vector<std::vector<int>> true_sources = { {0,0,0}, {-3,0,-2} , {-3,-3,0}, {-3,-1,-1} }; //the last two are not strictly needed, but produce shorter names
	std::vector<std::string> source_names = { "1" , "s3", "w3", "x1" }; //names of the sources above

	//Go from S^{-5*sigma-5*lambda} to S^{5*sigma+5*lambda}. The basic irreducibles are chosen together with their names
	Factorization<rank_t, diff_t> F({ -5,-5 }, { 5,5 }, { {0,1,0},{2,2,0},{0,0,1},{2,0,1} }, { "asigma", "u2sigma", "alambda", "ulambda" });
	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	std::vector<std::string> names;
	names.reserve(F.size);
	for (int i = 0; i < F.size; i++) {
		names.push_back(F.getname(i));
		std::cout << names[i] << "\n";
	}
	F.draw(names); //draws the graph
	return 0;
}
