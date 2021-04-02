#pragma once
#include <iostream>
#include <array>
#include <chrono>
#include "C4Verify.ipp"

namespace mackey
{
	namespace test
	{

		std::array<int, 4> interfaceAdd()
		{
			std::array<int, 4> ranges;
			bool flag = 1;
			while (flag)
			{
				for (int i = 0; i < 4; i++)
				{
					std::cin >> ranges[i];
				}
				if (ranges[0] > ranges[1] || ranges[2] > ranges[3])
				{
					std::cout << ranges[0] << "<= n <=" << ranges[1] << " and " << ranges[2] << "<= m<= " << ranges[3] << " is not a valid range. Please try again. \n";
					flag = 1;
				}
				else
					flag = 0;
			}
			std::cout << "Range " << ranges[0] << "<= n <=" << ranges[1] << " and " << ranges[2] << "<= m<= " << ranges[3] << " selected \n";
			std::cout << "Press any key to confirm and compute\n ";
			std::cin.get();
			return ranges;
		}

		std::array<int, 4> interfaceMul()
		{
			std::array<int, 4> ranges;
			bool flag = 1;
			while (flag)
			{
				for (int i = 0; i < 4; i++)
				{
					std::cin >> ranges[i];
				}
				if (ranges[0] < 0 || ranges[1] < 0 || ranges[2] < 0 || ranges[3] < 0)
				{
					std::cout << " All entered numbers must be nonnegative. Please try again. \n";
					flag = 1;
				}
				else
					flag = 0;
			}
			std::cout << "Max powers selected are asigma^" << ranges[0] << " , u2sigma^" << ranges[1] << " , alambda^" << ranges[2] << " and ulambda^" << ranges[3] << ".\n";
			std::cout << "Press any key to confirm and compute\n ";
			std::cin.get();
			return ranges;
		}

		template <typename group_t>
		void test_range(bool add, bool mul, bool silence, bool performance_test, const std::array<int, 4> &rangeAdd, const std::array<int, 4> &rangeMul)
		{

			if (add)
			{
				auto begin = std::chrono::high_resolution_clock::now();
				C4MackeyTest<group_t>(rangeAdd[0], rangeAdd[1], rangeAdd[2], rangeAdd[3], silence);
				auto end = std::chrono::high_resolution_clock::now();
				if (performance_test)
					std::cout << "Additive completed in: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "sec" << std::endl;
			}
			if (mul)
			{
				auto begin = std::chrono::high_resolution_clock::now();
				C4Multtest<group_t>(rangeMul[0], rangeMul[1], rangeMul[2], rangeMul[3], silence);
				auto end = std::chrono::high_resolution_clock::now();
				if (performance_test)
					std::cout << "Multiplicative completed in: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "sec" << std::endl;
			}
		}

		template <typename group_t>
		void basic_C4_add_mul_Z_coeff_test(bool userguided, bool print, bool lowranges, bool performance_test)
		{

			if (userguided)
			{
				std::cout << "The RO(C4) homology of a point in Z coefficients \n\n";
				std::cout << "Please enter \"a\" (no quotes) if you want the additive structure and \"m\" (no quotes) if you want the multiplicative one\n";
				bool flag = 1;
				std::string choice;

				while (flag)
				{
					std::cin >> choice;
					if (choice != "a" && choice != "m")
					{
						std::cout << "Please enter \"a\" (no quotes) if you want the additive structure and \"m\" (no quotes) if you want the multiplicative one\n";
					}
					else
					{
						flag = 0;
					}
				}

				if (choice == "a")
				{
					std::cout << "Additive structure selected.\nTo compute H_kS^{n*sigma+ m*lambda} for a<=n<=b and c<=m<=d enter the a b c d separated by spaces\nExample: -3 0 -5 10\n";
					auto ranges = interfaceAdd();
					std::cin.get();
					test_range<group_t>(1, 0, 0, 0, ranges, ranges);
				}
				if (choice == "m")
				{
					std::cout << "Multiplicative structure selected.\nWe will write all generators in H_kS^{n*sigma+ m*lambda} in terms of Euler classes asigma, u2sigma, orientation classes alambda, ulambda, the transfers w3 and x11 and finally the generator s3\n";
					std::cout << "Please enter the maximum powers the asigma, u2sigma, alambda, ulambda can take (either in numerators or denominators) as a b c d respectively, separated by spaces.\nExample 3 2 1 5\n";
					auto ranges = interfaceMul();
					std::cin.get();
					test_range<group_t>(0, 1, 0, 0, ranges, ranges);
				}

				std::cout << "\n\n Computation successful. Press any key to exit.";
				std::cin.get();
			}

			else
			{
				std::array<int, 4> rangeAdd, rangeMul;
				if (lowranges)
				{
					rangeAdd = {-3, 3, -3, 3};
					rangeMul = {3, 3, 3, 3};
				}
				else if (performance_test)
				{
					rangeAdd = {-30, 30, -30, 30};
					rangeMul = {13, 13, 13, 13};
				}
				else
				{
					rangeAdd = {-10, 10, -10, 10};
					rangeMul = {10, 10, 10, 10};
				}
				test_range<group_t>(1, 1, !print, performance_test, rangeAdd, rangeMul);
			}
		}
	}
}