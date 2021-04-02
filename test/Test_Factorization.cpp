#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen won't work otherwise.
#endif

#ifdef _OPENMP
#define MACKEY_USE_OPEN_MP
#endif

//Eigen directives
//#define EIGEN_USE_MKL_ALL
//#define EIGEN_NO_DEBUG //disable Eigen assertions
//#define NDEBUG
#include "Coefficients/Z_n.hpp"
#include "Groups/C2n.hpp"
#include "Factorization/Factorization.hpp"
#include <fstream>

typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::SparseMatrix<short> Zdiff_t;
typedef Eigen::SparseMatrix<Z2> Z2diff_t;
typedef Z_group_t = C4<rank_t,Zdiff_t>;
typedef Z2_group_t= C4<rank_t,Z2diff_t>;


using namespace mackey;

template<typename T>
void printgrapher(const T& F) {
	std::ofstream file;
	file.open("graph.dot");
	file << F.graph;
	file.close();
	file.open("spanningtree.dot");
	file << F.shortest_paths;
	file.close();
}

template<typename group_t>
void test_factorizationZ(int maximum, bool printgraph) {
	const int power = group_t::power;
	std::vector<std::string> names = { "sigma","lambda" };
	std::vector<std::string> names2 = { "2sigma","lambda" };

	//basic Irreducibles for all C_{2^n}. 
	std::vector<std::vector<int>> basicIrreducibles;
	std::vector<std::string> basicIrreducibles_names;
	basicIrreducibles.reserve(2 * power);
	basicIrreducibles_names.reserve(2 * power);
	for (int i = 0; i < power; i++) {
		std::vector<int> irr(power + 1);
		irr[i + 1] = 1;
		basicIrreducibles.push_back(irr);
		if (power==2)
			basicIrreducibles_names.push_back("a" + names[i]);
		else
			basicIrreducibles_names.push_back("a" + std::to_string(1 << (i + 1)));
	}
	for (int i = 0; i < power; i++) {
		std::vector<int> irr(power + 1);
		irr[i + 1] = 1;
		irr[0] = 2;
		if (i == 0)
			irr[1] = 2;
		basicIrreducibles.push_back(irr);
		if (power == 2)
			basicIrreducibles_names.push_back("u" + names2[i]);
		else
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

	Factorization<group_t> F(power, std::vector<int>(power, -maximum), std::vector<int>(power, maximum), basicIrreducibles, basicIrreducibles_names);
	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	//F.pass_disconnected(1);
	std::cout << F;
	if (printgraph)
		printgrapher(F);
}

template<typename group_t>
void test_factorizationZ2(int maximum, bool printgraph) {

	if (group_t::power != 2){
		std::cerr << "Only implemented for C4 currently";
		abort();
	}
	std::vector<std::vector<int>> true_sources = { {0,0,0}, {-3,0,-2} , {-2,0,-1} };
	std::vector<std::string> source_names = { "1" , "s", "theta" };
	Factorization<group_t> F(2, { -maximum,-maximum }, { maximum,maximum }, { {0,1,0},{1,1,0},{0,0,1},{2,0,1} }, { "asigma", "usigma", "alambda", "ulambda" });

	F.compute_with_sources(true_sources, source_names); //computes the factorizations
	std::cout << F;
	if (printgraph)
		printgrapher(F);
}

template<typename group_t>
void test_factorization(int maximum, bool printgraph1, bool printgraph2) {
	if constexpr (!SFINAE::is_finite_cyclic<typename group_t::diff_t::Scalar>::value)
		test_factorizationZ<group_t>(maximum, printgraph1);
	else
		test_factorizationZ2<group_t>(maximum, printgraph2);
}
int main()
{
	std::cout << "\n Factorization with Z coefficients\n";
	test_factorization<Z_group_t>(8, 0);

	std::cout << "\n Factorization with Z2 coefficients\n";
	test_factorization<Z_group_t>(8, 0);
	return 0;
}