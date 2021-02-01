//#define EIGEN_USE_MKL_ALL

#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen/cereal won't work otherwise.
#endif

#include <iostream>
#include "Cerealizer.h"
#include "C2n_Implementation.h"
#include "omp.h"

#include <chrono>

using namespace Mackey;

typedef Eigen::Matrix<char, 1, -1> rank_t;
typedef Eigen::SparseMatrix<char, 0, long> diff_t;
//typedef Eigen::Matrix<char,-1,-1> diff_t;

template<typename T>
bool boxnumber(const T& deg) {
	int poscount=0;
	int negcount=0;
	for (int i = 1; i < deg.size(); i++) {
		if (deg[i] > 0)
			poscount++;
		if (deg[i] < 0)
			negcount++;
	}
	if (poscount && negcount)
		return 1;
	else
		return 0;
}


template<typename T>
void doround(Factorization<rank_t, diff_t>& F, const T& round, bool do_indeter) {

	omp_lock_t lock;
	omp_init_lock(&lock);

	std::vector<char> done(round.size());
	std::cout << "Computing " << round.size() << " many Massey products \n";
#pragma omp parallel for num_threads(12) schedule(dynamic)
	for (int e = 0; e < round.size(); e++)
	{
		int i = round[e][0];
		int j = round[e][1];
		int k = round[e][2];
		//std::cout << F.getdegree(i) << " , " << F.getdegree(j) << " , " << F.getdegree(k) << "\n";
		auto a = ROMassey<rank_t, diff_t>(power, F.getdegree(i), F.getdegree(j), F.getdegree(k), find(F.getelement(i), 1), find(F.getelement(j), 1), find(F.getelement(k), 1), do_indeter);

		if (!do_indeter)
			a.noIndeterminacy = 1;
		done[e] = 1;
		if (a.exists && a.noIndeterminacy && a.normalBasis.size()==1 && a.normalBasis[0]!=0) {
			auto deg = F.getdegree(i) + F.getdegree(j) + F.getdegree(k);
			deg[0] += 1;
			rank_t el(1);
			el << 1;
			omp_set_lock(&lock);
			std::string name;
			if (F.getelementindex(deg, el)!=-1)
				name = F.getname(F.getelementindex(deg, el));
			std::cout << "<" << F.getname(i) << " , " << F.getname(j) << " , " << F.getname(k) << "> = " << (int)a.normalBasis[0] << " * " << name <<"\n";
			omp_unset_lock(&lock);
		}
	}

}



void doMassey(Factorization<rank_t, diff_t>& F) {

	std::vector<std::array<int, 3>> allrounds, firstround, secondround, thirdround;

		std::vector<int> acceptable;
		std::vector<std::vector<int>> degrees;
		for (int i = 0; i < F.number_of_generators; i++) {
			degrees.push_back(F.getdegree(i));
			if (!F.getname(i).empty() && boxnumber(F.getdegree(i)) == 0)
				acceptable.push_back(i);
		}
		allrounds.reserve(F.number_of_generators);
		bool flag = 0;
		for (short i = 0; i < acceptable.size(); i++) {
			if (flag == 1)
				break;
			for (short j = i; j < acceptable.size(); j++) {
				if (flag == 1)
					break;
				for (short k = j; k < acceptable.size(); k++) {
					if (flag == 1)
						break;
					if (degrees[acceptable[i]][0] + degrees[acceptable[j]][0] + degrees[acceptable[k]][0] != -4)
						continue;
					bool flag2 = 0;
					for (int u=1; u<power; u++)
						if (degrees[acceptable[i]][u] + degrees[acceptable[j]][u] + degrees[acceptable[k]][u] != 0) {
							flag2 = 1;
							break;
						}
					if (flag2 ==1)
						continue;
					if (degrees[acceptable[i]][power] + degrees[acceptable[j]][power] + degrees[acceptable[k]][power] != -2)
						continue;
					allrounds.push_back({ acceptable[i],acceptable[j],acceptable[k] });
					//if (allrounds.size() == 500)
					//	flag = 1;
					if (j!=k)
						allrounds.push_back({ acceptable[i],acceptable[k],acceptable[j] });
					if (i!=j)
						allrounds.push_back({ acceptable[j],acceptable[i],acceptable[k] });
					if (i != j && i!=k)
						allrounds.push_back({ acceptable[j],acceptable[k],acceptable[i] });
					if (i != k)
						allrounds.push_back({ acceptable[k],acceptable[i],acceptable[j] });
					if (i != k)
						allrounds.push_back({ acceptable[k],acceptable[j],acceptable[i] });
				}
			}
		}
	
		std::vector<char> selectfirst(allrounds.size());
	#pragma omp parallel for num_threads(12) schedule(dynamic)
		for (int e = 0; e < allrounds.size(); e++)
		{
			int i = allrounds[e][0];
			int j = allrounds[e][1];
			int k = allrounds[e][2];
			auto deg_l = F.getdegree(i) + F.getdegree(j);
			deg_l = Reindex(deg_l);
			deg_l[0] += 1;
			deg_l = invReindex(deg_l);
			auto deg_r = F.getdegree(j) + F.getdegree(k);
			deg_r = Reindex(deg_r);
			deg_r[0] += 1;
			deg_r = invReindex(deg_r);
	
			auto ind_l = ROHomology<rank_t, diff_t>(power, deg_l);
			auto ind_r = ROHomology<rank_t, diff_t>(power, deg_r);
			if (ind_l.size() == 0 && ind_r.size() == 0)
				selectfirst[e]=1;
		}
	
		firstround.reserve(allrounds.size());
		secondround.reserve(allrounds.size());
		for (int e = 0; e < allrounds.size(); e++) {
			if (selectfirst[e] == 1)
				firstround.push_back(allrounds[e]);
			else
				secondround.push_back(allrounds[e]);
		}
	//loader(firstround, "round1", "binary");
	doround(F, firstround, 0);
	doround(F, secondround,1);
	//doround(F, thirdround);
}



int main() {/*


	MultiplicationTable<rank_t, diff_t> M570, M589;
	loader(M570, "batch99_570", "binary");
	loader(M589, "batch99_589", "binary");

	bool a,b,c,d,e,s;
	a = (M570.NonZeroHomology == M589.NonZeroHomology);
	for (int i = 0; i < M570.degree.size(); i++)
		if (M570.degree[i] != M589.degree[i])
			std::cout << "H";
	b = (M570.degree == M589.degree);
	c = (M570.antidegree == M589.antidegree);
	d = (M570.index_product == M589.index_product);
	e = (M570.Greens == M589.Greens);
*/


	//std::vector<int> first = std::vector<int>{ -5, -3, -1 } ;
	//std::vector<int> second = std::vector<int>{ 4, 2,1};
	//std::vector<int> third = std::vector<int>{-3, 1, -2 };
	//auto u = ROMassey<rank_t, diff_t, std::vector<int>>(power, first, second, third, 1);

	////auto a = ROMassey<rank_t, diff_t, std::vector<int>>(power, { 0,-1,1}, { -2,0,-1}, { 0,-1,1 }, 1);


	//Eigen::Matrix<float, -1, -1> Out;
	//loader(Out, "problematic", "binary");
	//auto OUT=diagonalize<decltype(Out), decltype(Out), decltype(Out)>(Out, 1, 1);
	//std::cout << OUT.verify();

	//std::vector<int> first = std::vector<int>{ -3, -1, -1 } -std::vector<int>{6,2, 2} ;
	//std::vector<int> second = std::vector<int>{ 2, 2,0} +std::vector<int>{6, 2, 2};
	//std::vector<int> third = std::vector<int>{ -3, -1, -1 };
	//auto u = ROMassey<rank_t, diff_t, std::vector<int>>(power, first, second, third, 1);

//	
//	/*
//	std::vector<int> first = std::vector<int>{ -7, -1, -1, -1, -1 } -std::vector<int>{14, 0, 2, 3, 2} +std::vector<int>{2, 0, 0, 1, 0};
//	std::vector<int> second = std::vector<int>{ 10, 2, 2, 2, 0 } +std::vector<int>{14, 0, 2, 3, 2}-std::vector<int>{2, 0, 0, 1, 0};
//	std::vector<int> third = std::vector<int>{ -7, -1, -1, -1, -1 };
//
//
//	auto deg_l = first + second;
//	deg_l = Reindex(deg_l);
//	deg_l[0] += 1;
//	deg_l = invReindex(deg_l);
//	auto deg_r = second + third;
//	deg_r = Reindex(deg_r);
//	deg_r[0] += 1;
//	deg_r = invReindex(deg_r);
//
//	auto ind_l = ROHomology<rank_t, diff_t>(power, deg_l);
//	auto ind_r = ROHomology<rank_t, diff_t>(power, deg_r);
//	if (ind_l.size() == 0 && ind_r.size() == 0)
//		auto u = ROMassey<rank_t, diff_t, std::vector<int>>(4, first, second, third, 1);*/
//	//std::vector<int> xover(power + 1,-2);
//	//xover[0] = -4 * power + 3;
//	//xover[1] = -1;
//	//std::vector<int> xwith(power + 1);
//	//xwith[0] = -5;
//	//xwith[1] = -1;
//	//xwith.back() = -2;
//	//std::vector<int> prod(power + 1, 2);
//	//prod[0] = 4 * power - 2;
//	//auto u = ROMassey<rank_t, diff_t>(power, xover, xwith, prod, 1);
//
//
//
//	//std::vector<int> x(power + 1, -1);
//	//x[0] = -2*power+1;
//	//std::vector<int> xwith(power + 1);
//	//xwith[0] = -5;
//	//xwith[1] = -1;
//	//xwith.back() = -2;
//	//std::vector<int> prod(power + 1, 1);
//	//prod[0] = 2*power;
//	//prod[1] = 2;
//	//std::cout << x + xwith + prod << "\n";
//	//auto u=ROMassey<rank_t, diff_t>(power, x, xwith, prod, 1);
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


	int maximum = 5;

	Factorization<rank_t, diff_t> F(power, std::vector<int>(power, -maximum), std::vector<int>(power, maximum), basicIrreducibles, basicIrreducibles_names);
	//MultiplicationTable<rank_t, diff_t> M;
	//loader(M, "batch99", "binary");
	//Factorization<rank_t, diff_t> F(M, basicIrreducibles_names);

	auto nonS = true_sources;
	auto nonS_names = source_names;

	nonS.pop_back();
	nonS_names.pop_back();
	nonS.erase(nonS.begin() + 1);
	nonS_names.erase(nonS_names.begin() + 1);

	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	auto begin = std::chrono::high_resolution_clock::now();

	doMassey(F);



	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
	return 0;
}