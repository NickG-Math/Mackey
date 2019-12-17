//#define EIGEN_USE_MKL_ALL

#include <iostream>
#include "C4_Implementation.h"
#include "Mackey/Factorization.h"
#include <chrono>

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::Matrix<char, -1, -1> diff_t;

int main() {


std::vector<std::vector<int>> true_sources = { {0,0,0},  {-3,0,-2} , {-3,-3,0}, {-3,-1,-1} }; //the last two are not strictly needed, but produce shorter names
std::vector<std::string> source_names = { "1" , "s3", "w3", "x11"};


	auto begin = std::chrono::high_resolution_clock::now();

	//Go from S^{-5*sigma-5*lambda} to S^{5*sigma+5*lambda}. The basic irreducibles are chosen together with their names

	Factorization<rank_t, diff_t> F(2, { -5,-5 }, { 5,5 }, { {0,1,0},{2,2,0},{0,0,1},{2,0,1} }, { "asigma", "u2sigma", "alambda", "ulambda" });

	//Factorization<rank_t, diff_t> F(1, { -5,-5 }, { 5,5 }, {{1,1,0},{0,0,1},{2,0,1} }, {"usigma", "alambda", "ulambda" }); //for the middle level


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
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
	return 0;
}
