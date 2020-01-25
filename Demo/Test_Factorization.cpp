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
typedef Eigen::SparseMatrix<char, 0, long> diff_t;
//typedef Eigen::Matrix<char,-1,-1> diff_t;

int main() {

	//basic Irreducibles for all C_{2^n}. 
	std::vector<std::vector<int>> basicIrreducibles;
	std::vector<std::string> basicIrreducibles_names;
	basicIrreducibles.reserve(2 * power);
	basicIrreducibles_names.reserve(2 * power);
	for (int i = 0; i < power; i++) {
		std::vector<int> irr(power + 1);
		irr[i + 1] = 1;
		basicIrreducibles.push_back(irr);
		basicIrreducibles_names.push_back("a" + std::to_string(1 << (i + 1)));
	}
	for (int i = 0; i < power; i++) {
		std::vector<int> irr(power + 1);
		irr[i + 1] = 1;
		irr[0] = 2;
		if (i == 0)
			irr[1] = 2;
		basicIrreducibles.push_back(irr);
		basicIrreducibles_names.push_back("u" + std::to_string(1 << (i + 1)));
	}

	//sources for all C_{2^n}. 
	std::vector<std::vector<int>> true_sources(power + 2);
	true_sources[0].resize(power + 1);
	true_sources[1].resize(power + 1);
	true_sources[1][0] = -3;
	true_sources[1][1] = -3;
	for (int i = 2; i < power + 1; i++) {
		true_sources[i].resize(power + 1);
		for (int j = 1; j <= i; j++)
			true_sources[i][j] = -1;
		true_sources[i][0] = -2 * i + 1;
	}
	true_sources[power + 1].resize(power + 1);
	true_sources[power + 1][0] = -3;
	true_sources[power + 1][power] = -2;


	std::vector<std::string> source_names(power + 2);
	source_names[0] = "1";
	for (int i = 1; i < power + 1; i++) {
		source_names[i] = "x" + std::to_string(i);
	}
	source_names[power + 1] = "s";




	int maximum = 3;
	auto begin = std::chrono::high_resolution_clock::now();


	Factorization<rank_t, diff_t> F(power, std::vector<int>(power, -maximum), std::vector<int>(power, maximum), basicIrreducibles, basicIrreducibles_names);
	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	F.pass_disconnected(1);

	std::vector<std::string> names;
	std::vector<int> notfound;
	names.reserve(F.number_of_nodes);
	notfound.reserve(F.number_of_nodes);
	for (int i = 0; i < F.number_of_nodes; i++) {
		auto name = F.getname(i);
		if (name == "") {
			names.push_back("NOT FOUND");
			notfound.push_back(i);
		}
		else {
			names.push_back(name);
		}
		std::cout << i << " : " << names[i] << " at  " << F.getdegree(i)[0] << "," << F.getdegree(i)[1] << "," << F.getdegree(i)[2] << " and element: " << F.getelement(i) << "\n";
	}
	for (const auto& i : notfound) {
		std::cout << " Not found " << i<< " : "<< F.getdegree(i)[0] << "," << F.getdegree(i)[1] << "," << F.getdegree(i)[2] << " and element: " << F.getelement(i) << "\n";
	}

	F.draw(names); //draws the graph, make sure all nodes have names
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;

	return 0;
}
