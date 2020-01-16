//#define EIGEN_USE_MKL_ALL

#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen/cereal won't work otherwise.
#endif

#include <iostream>
#include "Mackey/Cerealizer.h"
#include "C4_Implementation.h"
#include "omp.h"

#include <chrono>

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::SparseMatrix<char, 0> diff_t;


void doMackey(Factorization<rank_t, diff_t>  F) {

	omp_lock_t lock;
	omp_init_lock(&lock);

	std::vector<std::array<int, 3>> tobedone;
	std::vector<std::vector<int>> degrees;

	for (int i = 0; i < F.number_of_generators; i++) {
		for (int j = 0; j < F.number_of_generators; j++) {
			for (int k = 0; k < F.number_of_generators; k++) {
				auto deg = F.getdegree(i) + F.getdegree(j) + F.getdegree(k);
				deg[0] += 1;
				if (deg[0] != -3 || deg[1] != 0 || deg[2] != -2)
					continue;
				tobedone.push_back({ i,j,k });
				degrees.push_back(deg);
			}
		}
	}
	std::cout << "Computing " << tobedone.size() << " many Massey products \n";
	rank_t one(1);
	one << 1;
#pragma omp parallel for num_threads(12) schedule(guided, 1)
	for (int e=0; e<tobedone.size(); e++)
	{
		int i = tobedone[e][0];
		int j = tobedone[e][1];
		int k = tobedone[e][2];
		auto deg = degrees[e];
		auto a = ROMassey<rank_t, diff_t, std::vector<int>>(2, F.getdegree(i), F.getdegree(j), F.getdegree(k), find(F.getelement(i), 1), find(F.getelement(j), 1), find(F.getelement(k), 1));

		omp_set_lock(&lock);

		if (a.exists && a.noIndeterminacy && a.basis.size() == 1 && a.basis[0] != 0) {
			std::cout << "[" << F.getname(i) << " , " << F.getname(j) << " , " << F.getname(k) << "] = " << (int)a.normalBasis[0] << " * " << F.getname(F.getelementindex(deg,one)) << "\n";
		}
		else if (a.exists && a.noIndeterminacy) {
			std::cout << "[" << F.getname(i) << " , " << F.getname(j) << " , " << F.getname(k) << "] = " << 0 << "\n";
		}
		omp_unset_lock(&lock);
	}

}


int main() {

	auto q = ROMassey< rank_t, Eigen::Matrix<Z<2>, -1, -1>, std::array<int,3>>(2, { 0,1,0 }, { -2,-2,0 }, { 1,1,0 });
	auto q1= ROMassey< rank_t, Eigen::Matrix<Z<2>, -1, -1>, std::array<int, 3>>(2, { 1,1,0 }, { -2,-2,0 }, { 0,1,0 });

	//std::vector<std::vector<int>> true_sources = { {0,0,0}, {1,0,1}, {-3,0,-2} , {-2,-2,0}, {-3,-1,-1}, {-2,0,-1} }; //the last two are not strictly needed, but produce shorter names
	//std::vector<std::string> source_names = { "1" , "q1", "s3", "w2", "x11", "v1" };

	std::vector<std::vector<int>> true_sources = { {0,0,0},  {-3,0,-2} , {-3,-3,0}, {-3,-1,-1} }; //the last two are not strictly needed, but produce shorter names
	std::vector<std::string> source_names = { "1" , "s3", "w3", "x11" };


	//Factorization<rank_t, diff_t> F(2, { -4,-4 }, { 4,4 }, { {0,1,0},{1,1,0},{0,0,1},{2,0,1} }, { "asigma", "usigma", "alambda", "ulambda" });

	Factorization<rank_t, diff_t> F(2, { -3,-3 }, { 3,3 }, { {0,1,0},{2,2,0},{0,0,1},{2,0,1} }, { "asigma", "u2sigma", "alambda", "ulambda" });

	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	auto begin = std::chrono::high_resolution_clock::now();

	doMackey(F);



	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
                                                    	return 0;
}