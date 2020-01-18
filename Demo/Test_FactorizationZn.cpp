//#define EIGEN_USE_MKL_ALL

#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen/cereal won't work otherwise.
#endif

#include <iostream>
#include "Mackey/Cerealizer.h"
#include "C2n_Implementation.h"
#include "Mackey/Factorization.h"
#include "Mackey/MultiplicationGraph_Connectivity.h"
#include <chrono>

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::SparseMatrix<char, 0> diff_t;

int main() {

	std::vector<std::vector<int>> basicIrreducibles;
	std::vector<std::string> basicIrreducibles_names;
	basicIrreducibles.reserve(2 * power);
	basicIrreducibles_names.reserve(2 * power);
	for (int i = 0; i < power; i++) {
		std::vector<int> irr(power+1);
		irr[i + 1] = 1;
		basicIrreducibles.push_back(irr);
		basicIrreducibles_names.push_back("a" + std::to_string(1 << (i+1)));
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

	std::vector<std::vector<int>> true_sources(power + 2);
	true_sources[0].resize(power+1);
	true_sources[1].resize(power + 1);
	true_sources[1][0] = -3;
	true_sources[1][1] = -3;
	for (int i = 2; i < power+1; i++) {
		true_sources[i].resize(power + 1);
		for (int j = 1; j <= i; j++)
			true_sources[i][j] = -1;
		true_sources[i][0] = -2 * i + 1;
	}
	true_sources[power+1].resize(power + 1);
	true_sources[power + 1][0] = -3;
	true_sources[power + 1][power] = -2;


	std::vector<std::string> source_names(power + 2);
	source_names[0] = "1";
	for (int i = 1; i < power+1; i++) {
		source_names[i]="x"+std::to_string(i);
	}
	source_names[power + 1] = "s";

	auto begin = std::chrono::high_resolution_clock::now();


	int maximum = 5;
	
	//MultiplicationTable<rank_t, diff_t> M(power, std::vector<int>(power, -maximum), std::vector<int>(power, maximum), basicIrreducibles);

	MultiplicationTable<rank_t, diff_t> M ;
	loader(M, "batch99", "binary");
	//M.extend({ -7,-7,-7,-7 }, { 7,7,7,7 });
//	saver(M,"Table16.bin","binary");
	//MultiplicationTable<rank_t, diff_t> M ;
	//loader(M, "Table.bin", "binary");
	//M.extend({ -10,-10,-10 }, { 10,10,10 });
	//saver(M, "Table.bin", "binary");



	MultiplicationGraphConnectivity<rank_t, diff_t> W(M);
	W.compute_with_sources(true_sources);
	Factorization<rank_t, diff_t> F(M, basicIrreducibles_names);
	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	//F.pass({ 10423,10424,10425 },1);
	//F.pass_unidentified(1);

	auto end = std::chrono::high_resolution_clock::now();

	std::vector<std::string> names;
	std::vector<int> notfound;
	names.reserve(F.number_of_generators);
	notfound.reserve(F.number_of_generators);
	for (int i = 0; i < F.number_of_generators; i++) {
		auto name = F.getname(i);
		if (name == "") {
			names.push_back("NOT FOUND");
			notfound.push_back(i);
		}
		else {
			names.push_back(name);
		}
		std::cout << i << " : " << names[i] << " at  ";
		auto u = F.getdegree(i);
		for (const auto & j:u )
			std::cout << j << ",";
		std::cout << " and element " << F.getelement(i) << "\n";
	}
	for (const auto &i :notfound) {
		auto u = F.getdegree(i);
		std::cout << i << " : " << " Not found ";
		for (const auto& j : u)
			std::cout << j << ",";
		std::cout << " and element " << F.getelement(i) << "\n";

	}
	F.compute_disconnected(F.number_of_generators);
	/////F.draw(names); //draws the graph
	/////F.draw(names,F.edgeid,basicIrreducibles_names);
	std::cout << "Time: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
	return 0;
}
