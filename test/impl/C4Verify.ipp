#pragma once
#include "Mackey_Functors/Additive.hpp"
#include <iostream>

///@file
///@brief Contains testing functions for the G=C4 case.
namespace mackey
{
	namespace test
	{

		//Human computed result to the RO(C_4) homology of a point over Z.
		template <typename group_t>
		std::vector<std::string> C4MackeyAnswer(int k, int n, int m)
		{
			static_assert(std::is_integral_v<typename group_t::diff_t::Scalar>, "C4MackeyAnswer only done for Z coefficients!");
			std::vector<std::string> MackeyName;
			std::vector<int> degree = {k, n, m};
			k = invReindex<group_t>(degree)[0];
			if (n >= 0 && m >= 0)
			{
				if (k == n + 2 * m && abs(n) % 2 == 0)
					MackeyName = {"Z", "111"};
				else if (k == n + 2 * m && abs(n) % 2 == 1)
					MackeyName = {"Z_-", "110"};
				else if (((abs(n) % 2 == 0 && k < n) || (abs(n) % 2 == 1 && k < n + 2 * m)) && abs(k) % 2 == 0 && k >= 0)
					MackeyName = {"Z/2", "002"};
				else if (abs(n) % 2 == 0 && k >= n && k < n + 2 * m && abs(k) % 2 == 0)
					MackeyName = {"Z/4", "024"};
				else if (abs(n) % 2 == 1 && k >= n && k < n + 2 * m && abs(k) % 2 == 1)
					MackeyName = {"overline Z/2", "020"};
				else
					MackeyName = {"0", "0"};
			}
			else if (n <= 0 && m <= 0)
			{
				if (k == n + 2 * m && n == -1 && m == 0)
					MackeyName = {"Z_-", "110"};
				else if (k == n + 2 * m && m == 0 && n != 0 && !(abs(n) % 2))
					MackeyName = {"p^*L", "111 # 1"};
				else if (k == n + 2 * m && n <= -3 && m == 0 && (abs(n) % 2) == 1)
					MackeyName = {"p^*L_-", "112 # 1"};
				else if (k == n + 2 * m && m != 0 && (abs(n) % 2) == 0)
					MackeyName = {"L", "111 # 01"};
				else if (k == n + 2 * m && m != 0 && (abs(n) % 2) == 1)
					MackeyName = {"L_-", "112 # 01"};
				else if (n < 0 && k <= 0 && abs(k) >= 3 && abs(k) % 2 == 1 && ((m == 0 && abs(k) < abs(n)) || (m != 0 && (abs(n) % 2) == 0 && abs(k) <= abs(n) + 1) || (m != 0 && (abs(n) % 2) == 1 && abs(k) < abs(n) + abs(2 * m))))
					MackeyName = {"Z/2", "002"};
				else if (abs(n) % 2 == 0 && k <= 0 && abs(k) >= abs(n) + 3 && abs(k) < abs(n + 2 * m) && abs(k) % 2 == 1)
					MackeyName = {"Z/4", "024"};
				else if (abs(n) % 2 == 1 && k <= 0 && abs(k) >= abs(n) + 3 && abs(k) < abs(n + 2 * m) && abs(k) % 2 == 0)
					MackeyName = {"overline Z/2", "020"};
				else
					MackeyName = {"0", "0"};
			}
			else if (n > 0 && m < 0)
			{
				m = -m;
				if (n % 2 == 1)
				{
					if ((n - 2 * m < k && k <= n - 4 && abs(k) % 2 == 1) || (0 <= k && k < n - 2 * m && abs(k) % 2 == 0))
						MackeyName = {"Z/2", "002"};
					else if (k == n - 2 * m && m >= 2)
						MackeyName = {"L_-", "112 # 01"};
					else if (k == n - 2 && m == 1)
						MackeyName = {"Z_-^flat", "110 # 0"};
					else if (k == n - 3 && n >= 3 && m >= 2)
						MackeyName = {"Q^sharp", "022"};
					else if ((n - 2 * m < k && k <= n - 5 && k < 0 && abs(k) % 2 == 0) || (k == -2 && n == 1 && m >= 2))
						MackeyName = {"overline Z/2", "020"};
					else if (n - 2 * m < k && k <= n - 5 && 0 <= k && abs(k) % 2 == 0)
						MackeyName = {"Z/2+overline Z/2", "020 + 002", "002 + 020"};
					else
						MackeyName = {"0", "0"};
				}
				else
				{
					if (0 <= k && k <= n - 4 && abs(k) % 2 == 0 && k != n - 2 * m)
						MackeyName = {"Z/2", "002"};
					else if (n - 2 * m < k && k < n - 3 && abs(k) % 2 == 1)
						MackeyName = {"Z/4", "024"};
					else if (k == n - 2 * m && n - 2 * m < 0 && m >= 2)
						MackeyName = {"L", "111 # 01"};
					else if (k == n - 2 && m == 1)
						MackeyName = {"L^sharp", "111 # 0"};
					else if (k == n - 3 && m >= 2)
						MackeyName = {"Q^sharp", "022"};
					else if (k == n - 2 * m && n - 2 * m >= 0 && m >= 2)
						MackeyName = {"Z/2+L", "002 + 111 # 01", " 111 # 01 + 002"};
					else
						MackeyName = {"0", "0"};
				}
			}
			else
			{
				n = -n;
				if (n % 2 == 0)
				{
					if (k == -n + 2 * m)
						MackeyName = {"Z", "111"};
					else if (-n + 1 <= k && k <= -3 && abs(k) % 2 == 1)
						MackeyName = {"Z/2", "002"};
					else if (-n + 2 <= k && k < -n + 2 * m && abs(k) % 2 == 0)
						MackeyName = {"Z/4", "024"};
					else if (k == -n)
						MackeyName = {"Q", "022 # 1"};
					else
						MackeyName = {"0", "0"};
				}
				else
				{
					if (k == -n + 2 * m && k >= -1)
						MackeyName = {"Z_-", "110"};
					else if ((-n + 1 <= k && k < 2 * m - n && abs(k) % 2 == 0) || (2 * m - n < k && k <= -3 && abs(k) % 2 == 1 && k != 2 * m - n))
						MackeyName = {"Z/2", "002"};
					else if (-1 <= k && k < -n + 2 * m && abs(k) % 2 == 1)
						MackeyName = {"overline Z/2", "020"};
					else if (k == -n && k <= -3)
						MackeyName = {"Q", "022 # 1"};
					else if (k >= -n + 2 && k < -n + 2 * m && k <= -3 && abs(k) % 2 == 1)
						MackeyName = {"Z/2+overline Z/2", "020 + 002", "002 + 020"};
					else if (k == -n + 2 * m && k <= -3)
						MackeyName = {"Z/2+Z_-", "002 + 110", "110 + 002"};
					else
						MackeyName = {"0", "0"};
				}
			}
			return MackeyName;
		}

		template <typename T>
		void print_to_stringstream(std::stringstream &ss, T v)
		{
			ss << v;
		}

		template <typename T, typename... Args>
		void print_to_stringstream(std::stringstream &ss, T v, Args... args)
		{
			ss << v;
			print_to_stringstream(ss, args...);
		}

		template <typename T, typename... Args>
		void print_Verification_or_Crash(bool silence, bool condition, const T &negative, const Args... affirmative)
		{
			std::stringstream ss;
			print_to_stringstream(ss, affirmative...);
			if (condition && !silence)
				std::cout << "Verified: " << ss.str() << "\n";
			else if (!condition)
			{
				std::cerr << "Not Verified! Should be " << ss.str() << "\n but instead got: " << negative;
				abort();
			}
		}


		///Tests the additive structure in the range rangeNLow<=n<=rangeNHigh and rangeMLow<=m<=rangeMHigh. silence=1 if we don't want console output.
		template <typename group_t>
		void C4MackeyTest(int rangeNLow, int rangeNHigh, int rangeMLow, int rangeMHigh, bool silence)
		{
			std::vector<int> minsphere = {rangeNLow, rangeMLow};
			std::vector<int> maxsphere = {rangeNHigh, rangeMHigh};
			AdditiveStructure<group_t> A(minsphere, maxsphere);
			if (!silence)
			{
				InterpolatingVectorGenerator<std::vector<int>> spheres(minsphere, maxsphere);
				for (const auto &i : spheres)
				{
					auto n = i[0];
					auto m = i[1];
					auto vector = A.getMackey(i);
					int k = 0;
					for (const auto &j : vector)
					{
						auto name = j.name;
						auto MackeyName = C4MackeyAnswer<group_t>(k, n, m);
						auto reindexed = invReindex<group_t>({k, n, m})[0];
						print_Verification_or_Crash(silence, MackeyName[1] == name || (MackeyName.size() == 3 && MackeyName[2] == name), name, "The k=", reindexed, " homology of the n=", n, " and m=", m, " sphere is ", MackeyName[0]);
						k++;
					}
				}
			}
		}

		///Tests the additive structure using powers of asigma, u2sigma, alambda, ulambda in a symmetric range described by the rangeN1,rangeN2,rangeM1,rangeM2 and rangeMLow<=m<=rangeMHigh respectively. silence=1 if we don't want console output.
		template <typename group_t>
		void C4Multtest(int rangeN1, int rangeN2, int rangeM1, int rangeM2, bool silence)
		{
			typedef typename group_t::rank_t rank_t;
			rank_t basis, element, elementnow, elementless;
			rank_t asigma(3), u2sigma(3), alambda(3), ulambda(3), w3(3), x11(3), s3(3);
			asigma << 0, 1, 0;
			u2sigma << 2, 2, 0;
			alambda << 0, 0, 1;
			ulambda << 2, 0, 1;
			w3 << -3, -3, 0;
			x11 << -3, -1, -1;
			s3 << -3, 0, -2;

			for (int n1 = 0; n1 <= rangeN1; n1++)
				for (int n2 = 0; n2 <= rangeN2; n2++)
					for (int m1 = 0; m1 <= rangeM1; m1++)
						for (int m2 = 0; m2 <= rangeM2; m2++)
							if (n1 <= 1 || m2 == 0)
							{
								elementnow = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda;
								if (n1 == 0 && n2 == 0 && m1 == 0 && m2 == 0)
									continue;
								if (n1 > 0)
								{
									elementless = elementnow - asigma;
									basis = ROGreen<group_t, rank_t>(2, asigma, elementless);
								}
								else if (n2 > 0)
								{
									elementless = elementnow - u2sigma;
									basis = ROGreen<group_t, rank_t>(2, u2sigma, elementless);
								}
								else if (m1 > 0)
								{
									elementless = elementnow - alambda;
									basis = ROGreen<group_t, rank_t>(2, alambda, elementless);
								}
								else if (m2 > 0)
								{
									elementless = elementnow - ulambda;
									basis = ROGreen<group_t, rank_t>(2, ulambda, elementless);
								}
								print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == 1, basis, "generator asigma^", n1, " * u2sigma^", n2, " * alambda^", m1, " * ulambda^", m2);
							}

			basis = ROGreen<group_t, rank_t>(2, 2 * asigma, ulambda);
			print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == 2, basis, "gold relation: ", "asigma^2 * ulambda = 2 * u2sigma * alambda");
			for (int n2 = 0; n2 >= -rangeN2; n2--)
				for (int m2 = 0; m2 >= -rangeM2; m2--)
				{
					element = n2 * u2sigma + m2 * ulambda;
					if (n2 == 0 && m2 == 0)
						continue;
					else if (n2 < 0)
						basis = ROGreen<group_t, rank_t>(2, u2sigma, element);
					else
						basis = ROGreen<group_t, rank_t>(2, ulambda, element);
					int k = 1;
					if (n2 == -1 && m2 == 0)
						k = 2;
					else if (n2 == 0 && m2 == -1)
						k = 4;
					int l = 4;
					if (m2 == 0)
						l = 2;
					print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == k, basis, "generator ", l, "/u2sigma^", -n2);
				}

			for (int n2 = 0; n2 >= -rangeN2; n2--)
				for (int m1 = 0; m1 >= -rangeM1; m1--)
					for (int m2 = 0; m2 >= -rangeM2; m2--)
					{
						element = n2 * u2sigma + m1 * alambda + m2 * ulambda + s3;
						if (n2 == 0 && m1 == 0 && m2 == 0)
							continue;
						if (n2 < 0)
							basis = ROGreen<group_t, rank_t>(2, u2sigma, element);
						else if (m1 < 0)
							basis = ROGreen<group_t, rank_t>(2, alambda, element);
						else if (m2 < 0)
							basis = ROGreen<group_t, rank_t>(2, ulambda, element);
						print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == 1, basis, "generator s3/(u2sigma^", -n2, "* alambda^", -m1, " * ulambda^", -m2, ")");
					}

			for (int n1 = 0; n1 >= -rangeN1; n1--)
				for (int n2 = 0; n2 >= -rangeN2; n2--)
					for (int m1 = 0; m1 >= -rangeM1; m1--)
						for (int m2 = 0; m2 >= -rangeM2; m2--)
						{
							if (n1 < 0 && m2 < 0)
								continue;
							if (n1 == 0 && n2 == 0 && m1 == 0 && m2 == 0)
								continue;
							element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda + x11;
							if (n1 < 0)
								basis = ROGreen<group_t, rank_t>(2, asigma, element);
							if (n2 < 0)
								basis = ROGreen<group_t, rank_t>(2, u2sigma, element);
							else if (m1 < 0)
								basis = ROGreen<group_t, rank_t>(2, alambda, element);
							else if (m2 < 0)
								basis = ROGreen<group_t, rank_t>(2, ulambda, element);
							print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == 1, basis, "generator x11/(asigma^", -n1, " * u2sigma^", -n2, " * alambda^", -m1, " * ulambda^", -m2, ")");
						}

			for (int n2 = -1; n2 >= -rangeN2; n2--)
				for (int m2 = 1; m2 <= rangeM2; m2++)
				{
					element = n2 * u2sigma + m2 * ulambda;
					basis = ROGreen<group_t, rank_t>(2, u2sigma, element);
					print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == 1, basis, "generator ulambda^", m2, " / u2sigma^", -n2);
				}

			for (int n1 = 0; n1 >= -1; n1--)
				for (int n2 = 0; n2 >= -rangeN2; n2--)
					for (int m1 = 1; m1 <= rangeM1; m1++)
						for (int m2 = 0; m2 <= rangeM2; m2++)
						{
							element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda;
							if (n1 < 0)
								basis = ROGreen<group_t, rank_t>(2, asigma, element);
							else if (n2 < 0)
								basis = ROGreen<group_t, rank_t>(2, u2sigma, element);
							int k = 1;
							if ((n1 == -1 && n2 == 0 && m2 == 0) || (n1 == -1 && m2 != 0) || (n1 == 0 && n2 == -1 && m2 == 0))
								k = 2;
							if (m2 == 0 || n1 != 0)
								print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == k, basis, "generator (2*alambda^", m1, " * ulambda^", m2, ")/(asigma^", -n1, " * u2sigma^", -n2, ")");
							else
								print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == k, basis, "generator (alambda^", m1, " * ulambda^", m2, ")/(asigma^", -n1, " * u2sigma^", -n2, ")");
						}

			for (int n1 = 0; n1 >= -rangeN1; n1--)
				for (int n2 = 0; n2 >= -rangeN2; n2--)
					for (int m1 = 1; m1 <= rangeM1; m1++)
					{
						element = n1 * asigma + n2 * u2sigma + m1 * alambda + w3;
						if (n1 == 0 && n2 == 0)
							basis = ROGreen<group_t, rank_t>(2, m1 * alambda, w3);
						if (n1 < 0)
							basis = ROGreen<group_t, rank_t>(2, asigma, element);
						else if (n2 < 0)
							basis = ROGreen<group_t, rank_t>(2, u2sigma, element);
						print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == 1, basis, "generator (alambda^", m1, " * w3)/(asigma^", -n1, " * u2sigma^", -n2, ")");
					}

			for (int n2 = 1; n2 <= rangeN2; n2++)
				for (int m2 = -1; m2 >= -rangeM2; m2--)
				{
					element = n2 * u2sigma + m2 * ulambda;
					basis = ROGreen<group_t, rank_t>(2, ulambda, element, 0, 1);
					if (m2 == -2 || m2 == -1)
					{
						bool condition1 = (basis.size() == 1 && basis(0) == 2) || (basis.size() == 2 && basis(0) == 0 && basis(1) == 1);
						bool condition2 = (basis.size() == 1 && basis(0) == 1) || (basis.size() == 2 && (basis(0) == 0 || basis(0) == 1) && basis(1) == 1);
						if (m2 == -2 || m2 == -1)
						{
							int k = 4;
							if (condition1 && m2 == -1)
								k = 2;
							print_Verification_or_Crash(silence, condition1 || condition2, basis, "generator ", k, "*u2sigma^", n2, " / ulambda^", -m2);
						}
					}
				}

			for (int n1 = 0; n1 <= 1; n1++)
				for (int n2 = 1; n2 <= rangeN2; n2++)
					for (int m1 = -1; m1 >= -rangeM1; m1--)
						for (int m2 = -1; m2 >= -rangeM2; m2--)
						{
							element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda + s3;
							if (m1 < 0)
								basis = ROGreen<group_t, rank_t>(2, alambda, element);
							else if (m2 < 0)
								basis = ROGreen<group_t, rank_t>(2, ulambda, element);
							print_Verification_or_Crash(silence, basis.size() == 1 && basis(0) == 1, basis, "generator (asigma^", n1, " * u2sigma^", n2, " * s3)/(alambda^", -m1, " * ulambda^", -m2, ")");
						}

			for (int n1 = 3; n1 <= rangeN1; n1++)
				for (int n2 = 1; n2 <= rangeN2; n2++)
					for (int m1 = -1; m1 >= -rangeM1; m1--)
					{
						element = n1 * asigma + n2 * u2sigma + m1 * alambda;
						basis = ROGreen<group_t, rank_t>(2, alambda, element);
						print_Verification_or_Crash(silence, (basis.size() == 1 && basis(0) == 1) || (basis.size() == 2 && basis(0) == 1 && basis(1) == 0), basis, "generator (asigma^", n1, " * u2sigma^", n2, ")/ alambda^", -m1);
					}
		}
	}
}