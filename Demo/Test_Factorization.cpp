//#define EIGEN_USE_MKL_ALL

#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen/cereal won't work otherwise.
#endif

#include <iostream>
#include "C4_Implementation.h"
#include "Mackey/Factorization.h"
#include <chrono>

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;

int main() {

	//typedef Eigen::Matrix<Z<2>, -1, -1> diff_t;
	//std::vector<std::vector<int>> true_sources = { {0,0,0}, {1,0,1}, {-3,0,-2} , {-2,-2,0}, {-3,-1,-1}, {-2,0,-1} };
	//std::vector<std::string> source_names = { "1" , "q1", "s3", "w2", "x11", "v1" };

	typedef Eigen::Matrix<char, -1, -1> diff_t;
	std::vector<std::vector<int>> true_sources = { {0,0,0},  {-3,-3,0}, {-3,-1,-1} }; //the last two are not strictly needed, but produce shorter names
	std::vector<std::string> source_names = { "1" , "w3", "x11" };


	auto begin = std::chrono::high_resolution_clock::now();

	//Go from S^{-5*sigma-5*lambda} to S^{5*sigma+5*lambda}. The basic irreducibles are chosen together with their names

	//Factorization<rank_t, diff_t> F(2, { -7,-7 }, { 7,7 }, { {0,1,0},{2,2,0},{0,0,1},{2,0,1} }, { "asigma", "u2sigma", "alambda", "ulambda" });

	//Factorization<rank_t, diff_t> F(2, { -4,-4 }, { 4,4 }, { {0,1,0},{1,1,0},{0,0,1},{2,0,1} }, { "asigma", "usigma", "alambda", "ulambda" });
	

	MultiplicationGraphConnectivity<rank_t, diff_t> F(2, { -7,-7 }, { 7,7 }, { {0,1,0},{2,2,0},{0,0,1},{2,0,1} });
	//MultiplicationGraphConnectivity<rank_t, diff_t> F(2, { -4,-4 }, { 4,4 }, { {0,1,0},{1,1,0},{0,0,1},{2,0,1} });


	//Factorization<rank_t, diff_t> F(1, { -5,-5 }, { 5,5 }, {{1,1,0},{0,0,1},{2,0,1} }, {"usigma", "alambda", "ulambda" }); //for the middle level


	F.compute_with_sources(true_sources); //computes the factorizations


	//F.compute_with_sources(true_sources, source_names); //computes the factorizations
	//F.pass_disconnected(1);

	//std::vector<std::string> names;
	//std::vector<int> notfound;
	//names.reserve(F.number_of_generators);
	//notfound.reserve(F.number_of_generators);
	//for (int i = 0; i < F.number_of_generators; i++) {
	//	auto name = F.getname(i);
	//	if (name == "") {
	//		names.push_back("NOT FOUND");
	//		notfound.push_back(i);
	//	}
	//	else {
	//		names.push_back(name);
	//	}
	//	std::cout << i << " : " << names[i] << " at  " << F.getdegree(i)[0] << "," << F.getdegree(i)[1] << "," << F.getdegree(i)[2] << " and element: " << F.getelement(i) << "\n";
	//}
	//for (const auto& i : notfound) {
	//	std::cout << " Not found " << i<< " : "<< F.getdegree(i)[0] << "," << F.getdegree(i)[1] << "," << F.getdegree(i)[2] << " and element: " << F.getelement(i) << "\n";
	//}
	//F.draw(names); //draws the graph, make sure all nodes have names
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;

	return 0;
}
