//#define EIGEN_USE_MKL_ALL

#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen/cereal won't work otherwise.
#endif

#include <iostream>
#include "Mackey/Cerealizer.h"
#include "C2n_Implementation.h"
#include "omp.h"

#include <chrono>

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::SparseMatrix<char, 0, long> diff_t;

template<typename T>
int boxnumber(const T& deg) {
	int bn = 0;
	for (int i = 1; i < deg.size() - 1; i++)
		if (deg[i] * deg[i + 1] < 0)
			bn++;
	return bn;
}






void doMackey(Factorization<rank_t, diff_t>  F) {

	omp_lock_t lock;
	omp_init_lock(&lock);

	std::vector<std::array<int, 3>> tobedone;
	std::vector<std::vector<int>> degrees;

	std::vector<int> s(power + 1);
	s[0] = -3;
	s.back() = -2;

	for (int i = 0; i < F.number_of_generators; i++) {
		if (boxnumber(F.getdegree(i)) > 0 || F.getname(i).empty())
			continue;
		for (int j = 0; j < F.number_of_generators; j++) {
			if (boxnumber(F.getdegree(j)) > 0 || F.getname(j).empty())
				continue;
			for (int k = 0; k < F.number_of_generators; k++) {
				if (boxnumber(F.getdegree(k)) > 0 || F.getname(k).empty())
					continue;
				auto deg = F.getdegree(i) + F.getdegree(j) + F.getdegree(k);
				deg[0] += 1;

				if (deg!= s)
					continue;
				tobedone.push_back({ i,j,k });
				degrees.push_back(deg);
			}
		}
	}
	std::cout << "Computing " << tobedone.size() << " many Massey products \n";
	rank_t one(1);
	one << 1;
#pragma omp parallel for num_threads(12) schedule(dynamic)
	for (int e=0; e<tobedone.size(); e++)
	{
		int i = tobedone[e][0];
		int j = tobedone[e][1];
		int k = tobedone[e][2];
		auto deg = degrees[e];
		auto a = ROMassey<rank_t, diff_t>(power, F.getdegree(i), F.getdegree(j), F.getdegree(k), find(F.getelement(i), 1), find(F.getelement(j), 1), find(F.getelement(k), 1));

		omp_set_lock(&lock);

		if (a.exists && a.noIndeterminacy && a.basis.size() == 1 && a.normalBasis[0] == 1) {
			std::cout << "[" << F.getname(i) << " , " << F.getname(j) << " , " << F.getname(k) << "] = " << (int)a.normalBasis[0] << " * " << F.getname(F.getelementindex(deg,one)) << "\n";
		}
		omp_unset_lock(&lock);
	}

}



int main() {



	//std::vector<int> xover(power + 1,-2);
	//xover[0] = -4 * power + 3;
	//xover[1] = -1;
	//std::vector<int> xwith(power + 1);
	//xwith[0] = -5;
	//xwith[1] = -1;
	//xwith.back() = -2;
	//std::vector<int> prod(power + 1, 2);
	//prod[0] = 4 * power - 2;
	//auto u = ROMassey<rank_t, diff_t>(power, xover, xwith, prod, 1);



	//std::vector<int> x(power + 1, -1);
	//x[0] = -2*power+1;
	//std::vector<int> xwith(power + 1);
	//xwith[0] = -5;
	//xwith[1] = -1;
	//xwith.back() = -2;
	//std::vector<int> prod(power + 1, 1);
	//prod[0] = 2*power;
	//prod[1] = 2;
	//std::cout << x + xwith + prod << "\n";
	//auto u=ROMassey<rank_t, diff_t>(power, x, xwith, prod, 1);
//

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


	int maximum = 2;

	Factorization<rank_t, diff_t> F(power, std::vector<int>(power, -maximum), std::vector<int>(power, maximum), basicIrreducibles, basicIrreducibles_names);

	auto nonS = true_sources;
	nonS.pop_back();
	auto nonS_names = source_names;
	nonS_names.pop_back();

	nonS.erase(nonS.begin() + 1);
	nonS_names.erase(nonS_names.begin()+1);

	F.compute_with_sources(nonS, nonS_names); //computes the factorizations
	auto begin = std::chrono::high_resolution_clock::now();

	doMackey(F);



	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
    return 0;
}