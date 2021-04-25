#pragma once
#include <iostream>
#include "Spaces/BC4S2.hpp"

namespace mackey
{
	namespace test
	{
		template <typename group_t>
		void test_BC4S2()
		{
			int dimension = 15;
			auto B2 = BC4S2<group_t>(dimension);
			auto bottomlevel = B2.ROHomology(0, { 0,0 });
			for (int i = 0; i < bottomlevel.size(); i++) {
				if (bottomlevel[i].number_of_summands() != 1) {
					std::cerr << "ERROR! ";
					abort();
				}
				else
					std::cout << "Verified bottom k^" << i << "(BC4S2) has dimension " << 1 << "\n";
			}
			std::cout << "\n";

			auto middlelevel = B2.ROHomology(1, { 0,0 });
			for (int i = 0; i < dimension/2; i++) {
				if (middlelevel[i].number_of_summands() != 2+i) {
					std::cerr << "ERROR! ";
					abort();
				}
				else
					std::cout << "Verified middle k^" << i << "(BC4S2) has dimension " << (2+i) << "\n";
			}
			std::cout << "\n";

			auto toplevel = B2.ROHomology(2, { 0,0 });
			for (int i = 0; i < dimension / 4; i++) {
				if (toplevel[i].number_of_summands() != 2*i+3) {
					std::cerr << "ERROR! ";
					abort();
				}
				else
					std::cout << "Verified top k^" << i << "(BC4S2) has dimension " << (2*i+3) << "\n";
			}
			std::cout << "\n";


			auto toplevel2sigma = B2.ROCohomology(2, { 2,0 });
			if (toplevel2sigma[3].number_of_summands() != 0) {
				std::cerr << "ERROR! ";
				abort();
			}
			else
				std::cout << "top k^{2sigma+1}=0 which verifies the nontrivial d2" << "\n";
		}

	}
}