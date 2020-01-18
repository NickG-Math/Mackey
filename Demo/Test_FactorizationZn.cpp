//#define EIGEN_USE_MKL_ALL

#ifdef _MSC_VER
#define _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING //Eigen won't work otherwise.
#endif

#include <iostream>
#include "C2n_Implementation.h"
#include "Mackey/Factorization.h"
#include <chrono>

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::Matrix<Z<2>, -1, -1> diff_t;

int main() {
	

	std::vector<std::vector<int>> true_sources = { {0,0,0}, {1,0,1}, {-3,0,-2} , {-2,-2,0}, {-3,-1,-1}, {-2,0,-1} }; 
	std::vector<std::string> source_names = { "1" , "q1", "s3", "w2", "x11", "v1" };


	auto begin = std::chrono::high_resolution_clock::now();


	//Go from S^{-5*sigma-5*lambda} to S^{5*sigma+5*lambda}. The basic irreducibles are chosen together with their names
	Factorization<rank_t, diff_t> F(2,{ -5,-5 }, {5,5}, { {0,1,0},{1,1,0},{0,0,1},{2,0,1} }, { "asigma", "usigma", "alambda", "ulambda" });
	//Factorization<rank_t, diff_t> F(1, { -5,-5 }, { 5,5 }, {{1,1,0},{0,0,1},{2,0,1} }, {"usigma", "alambda", "ulambda" }); //middle level


	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	std::vector<std::string> names;
	std::vector<int> notfound;
	names.reserve(F.size);
	notfound.reserve(F.size);
	for (int i = 0; i < F.size; i++) {
		auto name = F.getname(i);
		if (name == "") {
			names.push_back("NOT FOUND");
			notfound.push_back(i);
		}
		else {
			names.push_back(name);
		}
		std::cout << names[i] << " at  " << F.getdegree(i)[0] << "," << F.getdegree(i)[1] << "," << F.getdegree(i)[2]  << "\n";
	}
	for (const auto &i :notfound) {
		std::cout << " Not found " << F.getdegree(i)[0] << "," << F.getdegree(i)[1] << "," << F.getdegree(i)[2] << " and position " << F.getposition(i) << "\n";
	}
	F.draw(names); //draws the graph

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
	return 0;
}
