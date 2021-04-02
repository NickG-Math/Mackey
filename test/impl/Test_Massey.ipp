#pragma once
#include <iostream>
#include "Utility/OpenMP_Macros.hpp"
#include "Factorization/Factorization.hpp"
#include "Test_Factorization.ipp"

namespace mackey
{
	namespace test
	{

		template <typename T>
		bool involves_box_product(const T& deg)
		{
			int poscount = 0;
			int negcount = 0;
			for (auto d : deg)
			{
				poscount += (d > 0);
				negcount += (d < 0);
				if (poscount && negcount)
					return 1;
			}
			return 0;
		}

		template<typename group_t>
		std::vector<std::array<int, 3>> acceptable_triples(const Factorization<group_t>& F) {
			std::vector<std::array<int, 3>> triples;
			triples.reserve(F.number_of_generators * 3);
			for (int i = 0; i < F.number_of_generators; i++) {
				if (involves_box_product(F.getdegree(i)))
					continue;
				for (int j = 0; j < F.number_of_generators; j++) {
					if (involves_box_product(F.getdegree(j)))
						continue;
					for (int k = 0; k < F.number_of_generators; k++) {
						if (involves_box_product(F.getdegree(k)))
							continue;
						if (!F.getname(i).empty() && !F.getname(j).empty() && !F.getname(k).empty()) {
							auto deg = F.getdegree(i) + F.getdegree(j) + F.getdegree(k);
							deg[0]++;
							if (F.degreewithinrange(deg) && F.getname(deg).size() == 1 && !F.getname(deg).front().empty())
								triples.push_back({ i,j,k });
						}
					}
				}
			}
			return triples;
		}

		template<typename T>
		std::vector<T> operator*(T a, const std::vector<T>& v) {
			std::vector<T> y = v;
			for (auto& i : y)
				i *= a;
			return y;
		}

		template <typename group_t>
		void test_Massey() {
			auto F = getfactorizationZ<group_t>(4);
			auto triples = acceptable_triples(F);
			std::cout << "Computing " << triples.size() << " many Massey products\n ";
			int counter = 0;
			int counteraccess = 1;
			std::ofstream file;
			file.open("massey.txt");
			MACKEY_RUN_LOOP_PARALLEL
				for (int e = 0; e < triples.size(); e++)
				{
#pragma omp atomic
					counter++;
					if (counter == 1000)
					{
						MACKEY_RUN_BLOCK_SERIAL{
						if (counter == 1000) {
						std::cout << "\n\nDone " << (counteraccess * 1000) << " many Massey products\n\n";
						counter = 0;
						counteraccess++;
						}
						}
					}
					int i = triples[e][0];
					int j = triples[e][1];
					int k = triples[e][2];
					auto a = ROMassey<group_t>(group_t::power, F.getdegree(i), F.getdegree(j), F.getdegree(k));
					if (a.exists && a.noIndeterminacy && a.normalBasis.size() == 1 && a.normalBasis[0] != 0)
					{
						auto deg = F.getdegree(i) + F.getdegree(j) + F.getdegree(k);
						deg[0]++;
						auto name = F.getname(deg);
						std::stringstream ss;
						ss << "<" << F.getname(i) << " , " << F.getname(j) << " , " << F.getname(k) << "> = " << (int)a.normalBasis[0] << " * " << name << "\n";
						std::cout << ss.str();
						MACKEY_RUN_BLOCK_SERIAL
						{
							file << ss.str();
						}
					}
				}
			file.close();
		}
	}
}