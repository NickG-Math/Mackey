#ifdef _MSC_VER
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS //Eigen won't work otherwise.
#endif
#define MACKEY_USE_OPEN_MP
#include "Mackey.hpp"
#include <fstream>
// #include "Serialization/Cerealizer.hpp"

int main() {
	using namespace mackey;

	typedef C4<int64_t> group_t;

	AdditiveStructure<group_t> A({ -3, -4 }, { 5, 6 });
	std::cout << A << "\n";

	std::cout << "The identified Mackey functors are: " << A.identified() << "\n\n\n\n";

	std::vector<std::vector<int>> basic_irr = { {0, 1, 0}, {2, 2, 0}, {0, 0, 1}, {2, 0, 1} };
	std::vector<std::string> basic_irr_names = { "asigma", "u2sigma", "alambda", "ulambda" };

	auto F = Factorization<group_t>(2, { -5, -5 }, { 5, 5 }, basic_irr, basic_irr_names);

	F.compute_with_sources({ {0, 0, 0} }, { "1" });

	std::cout << "The generator(s) at degree 3,1,0 is (are): " << F.getname({ -2,-2,0 }) << "\n\n";

	std::cout << F << "\n\n";

	F.compute_with_sources({ {0, 0, 0}, {-3, 0, -2} }, { "1", "s" });

	std::cout << F << "\n\n";

	std::ofstream file;
	file.open("multgraph.dot");
	file << F.graph;
	file.close();

	file.open("shortestpaths.dot");
	file << F.shortest_paths;
	file.close();

	auto linear_combination = ROGreen<group_t>(2, { -4, -2, -1 }, { 2, 0,  1 });
	if (linear_combination.size() == 1)
		std::cout << F.getname({ -4,-2,-1 }) << " * " << F.getname({ 2,0,1 }) << " = " << linear_combination[0] << " * " << F.getname({ -2,-2,0 }) << "\n";

	auto Mass = ROMassey<group_t>(2, { 0, 0, 3 }, { -3, 0, -2 }, { 2, 0, 1 });
	if (Mass.noIndeterminacy && Mass.basis.size() == 1)
		std::cout << "<" << F.getname({ 0, 0, 3 }) << " , " << F.getname({ -3, 0, -2 }) << " , " << F.getname({ 2, 0, 1 }) << "> = " << Mass.basis[0] << " * " << F.getname({ 0,0,2 });


	typedef C8<int64_t> group2_t;


	AdditiveStructure<group2_t> A2({ -3, -3,-3 }, { 3, 3,3 });
	std::cout << A2 << "\n\n";
	std::cout << "The unknown 5 is =\n" << A2.unknown()[5].print() << "\n\n";


	//basic Irreducibles for C_8
	std::vector<std::vector<int>> basic_irr2 = { {0,1,0,0},{0,0,1,0},{0,0,0,1},{2,2,0,0},{2,0,1,0},{2,0,0,1} };
	std::vector<std::string> basic_names2 = { "a2","a4","a8","u2","u4","u8" };
	std::vector<std::vector<int>> sources2 = { {0,0,0,0}, {-3,0,0,-2} };
	std::vector<std::string> source_names2 = { "1", "s" };

	auto F2 = Factorization<group2_t>(3, std::vector{ -5,-5,-5 }, std::vector{ 5, 5, 5 }, basic_irr2, basic_names2);
	F2.compute_with_sources(sources2, source_names2);
	std::cout << F2.disconnected_degrees().size() << "\n\n"
		;
	F2.pass_unidentified();
	std::cout << F2.disconnected_degrees().size() << "\n\n";

	auto MC = MultConnectivity<group2_t>(F2);
	MC.compute_with_sources(sources2);
	std::cout << MC.disconnected_degrees.size() << "\n";


	// save(F,"filename");
	// MultTableData<group_t> M;
	// load(M,"filename");
	// Factorization<group_t> Fnew(std::move(M),basic_irr_names);
	return 0;
}

