<<<<<<< HEAD
#pragma once
#include <Mackey/Compute.h>
#include <string>
#include <iostream>

///@file
///@brief Contains testing functions for the G=C4 case.
namespace {

	std::string C4MackeyAnswer(int k, int n, int m) {
		std::string MackeyName;
		std::vector<int> sphere = { n,m };
		k = Mackey::invReindex(k, sphere);
		if (n >= 0 && m >= 0) {
			if (k == n + 2 * m && abs(n) % 2 == 0)
				MackeyName = "Z";
			else if (k == n + 2 * m && abs(n) % 2 == 1)
				MackeyName = "Z_-";
			else if (((abs(n) % 2 == 0 && k < n) || (abs(n) % 2 == 1 && k < n + 2 * m)) && abs(k) % 2 == 0 && k >= 0)
				MackeyName = "Z/2";
			else if (abs(n) % 2 == 0 && k >= n && k < n + 2 * m && abs(k) % 2 == 0)
				MackeyName = "Z/4";
			else if (abs(n) % 2 == 1 && k >= n && k < n + 2 * m && abs(k) % 2 == 1)
				MackeyName = "overline Z/2";
			else
				MackeyName = "0";
		}
		else if (n <= 0 && m <= 0) {
			if (k == n + 2 * m && n == -1 && m == 0)
				MackeyName = "Z_-";
			else if (k == n + 2 * m && m == 0 && n != 0 && !(abs(n) % 2))
				MackeyName = "p^*L";
			else if (k == n + 2 * m && n <= -3 && m == 0 && (abs(n) % 2) == 1)
				MackeyName = "p^*L_-";
			else if (k == n + 2 * m && m != 0 && (abs(n) % 2) == 0)
				MackeyName = "L";
			else if (k == n + 2 * m && m != 0 && (abs(n) % 2) == 1)
				MackeyName = "L_-";
			else if (n < 0 && k <= 0 && abs(k) >= 3 && abs(k) % 2 == 1 && ((m == 0 && abs(k) < abs(n)) || (m != 0 && (abs(n) % 2) == 0 && abs(k) <= abs(n) + 1) || (m != 0 && (abs(n) % 2) == 1 && abs(k) < abs(n) + abs(2 * m))))
				MackeyName = "Z/2";
			else if (abs(n) % 2 == 0 && k <= 0 && abs(k) >= abs(n) + 3 && abs(k) < abs(n + 2 * m) && abs(k) % 2 == 1)
				MackeyName = "Z/4";
			else if (abs(n) % 2 == 1 && k <= 0 && abs(k) >= abs(n) + 3 && abs(k) < abs(n + 2 * m) && abs(k) % 2 == 0)
				MackeyName = "overline Z/2";
			else
				MackeyName = "0";
		}
		else if (n > 0 && m < 0) {
			m = -m;
			if (n % 2 == 1) {
				if ((n - 2 * m < k && k <= n - 4 && abs(k) % 2 == 1) || (0 <= k && k < n - 2 * m && abs(k) % 2 == 0))
					MackeyName = "Z/2";
				else if (k == n - 2 * m && m >= 2)
					MackeyName = "L_-";
				else if (k == n - 2 && m == 1)
					MackeyName = "Z_-^flat";
				else if (k == n - 3 && n >= 3 && m >= 2)
					MackeyName = "Q^sharp";
				else if ((n - 2 * m < k && k <= n - 5 && k < 0 && abs(k) % 2 == 0) || (k == -2 && n == 1 && m >= 2))
					MackeyName = "overline Z/2";
				else if (n - 2 * m < k && k <= n - 5 && 0 <= k && abs(k) % 2 == 0)
					MackeyName = "Z/2+overline Z/2";
				else
					MackeyName = "0";
			}
			else {
				if (0 <= k && k <= n - 4 && abs(k) % 2 == 0 && k != n - 2 * m)
					MackeyName = "Z/2";
				else if (n - 2 * m < k && k < n - 3 && abs(k) % 2 == 1)
					MackeyName = "Z/4";
				else if (k == n - 2 * m && n - 2 * m < 0 && m >= 2)
					MackeyName = "L";
				else if (k == n - 2 && m == 1)
					MackeyName = "L^sharp";
				else if (k == n - 3 && m >= 2)
					MackeyName = "Q^sharp";
				else if (k == n - 2 * m && n - 2 * m >= 0 && m >= 2)
					MackeyName = "Z/2+L";
				else
					MackeyName = "0";
			}
		}
		else {
			n = -n;
			if (n % 2 == 0) {
				if (k == -n + 2 * m)
					MackeyName = "Z";
				else if (-n + 1 <= k && k <= -3 && abs(k) % 2 == 1)
					MackeyName = "Z/2";
				else if (-n + 2 <= k && k < -n + 2 * m && abs(k) % 2 == 0)
					MackeyName = "Z/4";
				else if (k == -n)
					MackeyName = "Q";
				else
					MackeyName = "0";
			}
			else {
				if (k == -n + 2 * m && k >= -1)
					MackeyName = "Z_-";
				else if ((-n + 1 <= k && k < 2 * m - n && abs(k) % 2 == 0) || (2 * m - n < k && k <= -3 && abs(k) % 2 == 1 && k != 2 * m - n))
					MackeyName = "Z/2";
				else if (-1 <= k && k < -n + 2 * m && abs(k) % 2 == 1)
					MackeyName = "overline Z/2";

				else if (k == -n && k <= -3)
					MackeyName = "Q";
				else if (k >= -n + 2 && k < -n + 2 * m && k <= -3 && abs(k) % 2 == 1)
					MackeyName = "Z/2+overline Z/2";

				else if (k == -n + 2 * m && k <= -3)
					MackeyName = "Z/2+Z_-";
				else
					MackeyName = "0";
			}
		}
		return MackeyName;
	}

}

///The namespace for testing the answers produced by Mackey (in case of G=C4).
namespace C4Test{

	using namespace Mackey;

	///Tests the additive structure in the range rangeNLow<=n<=rangeNHigh and rangeMLow<=m<=rangeMHigh. silence=1 if we don't want console output.
	template<typename rank_t, typename diff_t>
	void C4MackeyTest(int rangeNLow, int rangeNHigh, int rangeMLow, int rangeMHigh, bool silence)
	{
		std::vector<int> minimum = { rangeNLow, rangeMLow };
		std::vector<int> maximum = { rangeNHigh, rangeMHigh };
		auto degrees = DegreeConstruction(minimum, maximum);
		for (const auto& i : degrees)
		{
			auto n = i[0];
			auto m = i[1];
			auto M = ROHomology<rank_t, diff_t>(i);
			int k = 0;
			for (auto& j : M) {
				if (!silence) {
					auto name = identify_Mackey<rank_t, diff_t>(j);
					bool correct=1;
					if constexpr std::is_integral_v<typename diff_t::Scalar>{
						auto MackeyName = C4MackeyAnswer(k, n, m);
						if (MackeyName != name)
						throw(0);
					}
					std::cout << "The k=" << invReindex<std::vector<int>>(k, { n,m }) << " homology of the n=" << n << " and m=" << m << " sphere is " << name << "\n";
				}
				k++;
			}
		}
	}



	///Tests the additive structure using powers of asigma, u2sigma, alambda, ulambda in a symmetric range described by the rangeN1,rangeN2,rangeM1,rangeM2 and rangeMLow<=m<=rangeMHigh respectively. silence=1 if we don't want console output.
	template<typename rank_t, typename diff_t>
	void C4Multtest(int rangeN1, int rangeN2, int rangeM1, int rangeM2, bool silence) {
		rank_t basis, element, elementnow, elementless;
		rank_t asigma(3), u2sigma(3), alambda(3), ulambda(3), w3(3), x11(3), s3(3);
		asigma << 0, 1, 0; u2sigma << 2, 2, 0; alambda << 0, 0, 1; ulambda << 2, 0, 1; w3 << -3, -3, 0; x11 << -3, -1, -1; s3 << -3, 0, -2;


		for (int n1 = 0; n1 <= rangeN1; n1++) {
			for (int n2 = 0; n2 <= rangeN2; n2++) {
				for (int m1 = 0; m1 <= rangeM1; m1++) {
					for (int m2 = 0; m2 <= rangeM2; m2++) {
						if (n1 <= 1 || m2 == 0) {
							elementnow = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda;
							if (n1 == 0 && n2 == 0 && m1 == 0 && m2 == 0) {
								continue;
							}
							if (n1 > 0) {
								elementless = elementnow - asigma;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, elementless);
							}
							else if (n2 > 0) {
								elementless = elementnow - u2sigma;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, elementless);
							}
							else if (m1 > 0) {
								elementless = elementnow - alambda;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, elementless);
							}
							else if (m2 > 0) {
								elementless = elementnow - ulambda;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, elementless);
							}
							if (!silence) {
								if (basis.size() == 1 && basis(0) == 1) {
									std::cout << "Top Generator verified asigma^" << n1 << "*u2sigma^" << n2 << "*alambda^" << m1 << "*ulambda^" << m2 << "\n";
								}
								else {
									throw(1);
								}
							}
						}
					}
				}
			}
		}


		basis = ROGreen<rank_t, diff_t, rank_t>(2, 2 * asigma, ulambda);
		if (!silence) {
			if (basis.size() == 1 && basis(0) == 2) {
				std::cout << "Gold Verified!" << "\n";

			}
			else {
				throw(2);
			}
		}

		//element = -8 * u2sigma +6 * ulambda;
		//basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);



		for (int n2 = 0; n2 >= -rangeN2; n2--) {
			for (int m2 = 0; m2 >= -rangeM2; m2--) {
				element = n2 * u2sigma + m2 * ulambda;
				if (n2 == 0 && m2 == 0) {
					continue;
				}
				else if (n2 < 0) {
					basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
				}
				else {
					basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
				}
				if (!silence) {

					if (n2 == -1 && m2 == 0) {
						if (basis.size() == 1 && basis(0) == 2) {
							std::cout << "Top Generator verified: 2 / u2sigma^" << (-n2) << "\n";
						}
						else {
							throw(3);
						}
					}
					else if (n2 == 0 && m2 == -1) {
						if (basis.size() == 1 && basis(0) == 4) {
							std::cout << "Top Generator verified: 4 / ulambda^" << (-m2) << "\n";
						}
						else {
							throw(3);
						}
					}
					else {
						if (basis.size() == 1 && basis(0) == 1) {
							if (m2 == 0) {
								std::cout << "Top Generator verified: 2 / u2sigma^" << (-n2) << "\n";
							}
							else {
								std::cout << "Top Generator verified: 4 / u2sigma^" << (-n2) << "*ulambda^" << (-m2) << "\n";
							}
						}
						else {
							throw(3);
						}
					}
				}
			}
		}

		for (int n2 = 0; n2 >= -rangeN2; n2--) {
			for (int m1 = 0; m1 >= -rangeM1; m1--) {
				for (int m2 = 0; m2 >= -rangeM2; m2--) {
					element = n2 * u2sigma + m1 * alambda + m2 * ulambda + s3;
					if (n2 == 0 && m1 == 0 && m2 == 0) {
						continue;
					}
					if (n2 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
					}
					else if (m1 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
					}
					else if (m2 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
					}
					if (!silence) {

						if (basis.size() == 1 && basis(0) == 1) {
							std::cout << "Top Generator verified s3 / u2sigma^" << (-n2) << "*alambda^" << (-m1) << "*ulambda^" << (-m2) << "\n";
						}
						else {
							throw(4);
						}
					}
				}
			}
		}

		for (int n1 = 0; n1 >= -rangeN1; n1--) {
			for (int n2 = 0; n2 >= -rangeN2; n2--) {
				for (int m1 = 0; m1 >= -rangeM1; m1--) {
					for (int m2 = 0; m2 >= -rangeM2; m2--) {
						if (n1 < 0 && m2 < 0) {
							continue;
						}
						if (n1 == 0 && n2 == 0 && m1 == 0 && m2 == 0) {
							continue;
						}
						element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda + x11;
						if (n1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, element);
						}
						if (n2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
						}
						else if (m1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
						}
						else if (m2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
						}
						if (!silence) {

							if (basis.size() == 1 && basis(0) == 1) {
								std::cout << "Top Generator verified x11 / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "*alambda^" << (-m1) << "*ulambda^" << (-m2) << "\n";
							}
							else {
								throw(5);
							}
						}
					}
				}
			}
		}

		for (int n2 = -1; n2 >= -rangeN2; n2--) {
			for (int m2 = 1; m2 <= rangeM2; m2++) {
				element = n2 * u2sigma + m2 * ulambda;
				basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
				if (!silence) {

					if (basis.size() == 1 && basis(0) == 1) {
						std::cout << "Top Generator verified ulambda^" << (m2) << " / u2sigma^" << (-n2) << "\n";
					}
					else {
						throw(1);
					}
				}
			}
		}



		for (int n1 = 0; n1 >= -1; n1--) {
			for (int n2 = 0; n2 >= -rangeN2; n2--) {
				for (int m1 = 1; m1 <= rangeM1; m1++) {
					for (int m2 = 0; m2 <= rangeM2; m2++) {
						element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda;
						if (n1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, element);
						}
						else if (n2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
						}
						if (!silence) {

							if ((n1 == -1 && n2 == 0 && m2 == 0) || (n1 == -1 && m2 != 0) || (n1 == 0 && n2 == -1 && m2 == 0)) {
								if (basis.size() == 1 && basis(0) == 2) {
									if (m2 == 0 || n1 != 0) {
										std::cout << "Top Generator verified 2*alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";
									}
									else {
										std::cout << "Top Generator verified alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";
									}
								}
								else {
									throw(6);
								}
							}
							else {
								if (basis.size() == 1 && basis(0) == 1) {
									if (m2 == 0 || n1 != 0) {
										std::cout << "Top Generator verified 2*alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << " \n";
									}
									else {
										std::cout << "Top Generator verified alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";

									}
								}
								else {
									throw(6);
								}
							}
						}
					}
				}
			}
		}



		for (int n1 = 0; n1 >= -rangeN1; n1--) {
			for (int n2 = 0; n2 >= -rangeN2; n2--) {
				for (int m1 = 1; m1 <= rangeM1; m1++) {
					element = n1 * asigma + n2 * u2sigma + m1 * alambda + w3;
					if (n1 == 0 && n2 == 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, m1 * alambda, w3);
					}
					if (n1 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, element);
					}
					else if (n2 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
					}
					if (!silence) {

						if (basis.size() == 1 && basis(0) == 1) {
							std::cout << "Top Generator verified alambda^" << m1 << "*w3" << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";

						}
						else {
							throw(7);
						}
					}
				}
			}
		}

//element = 9 * u2sigma -10 * ulambda;
//basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element, 0, 1);


		for (int n2 = 1; n2 <= rangeN2; n2++) {
			for (int m2 = -1; m2 >= -rangeM2; m2--) {
				element = n2 * u2sigma + m2 * ulambda;
				basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element, 0, 1);
				if (!silence) {

					if (m2 == -2 || m2 == -1) {
						if ((basis.size() == 1 && basis(0) == 2) || (basis.size() == 2 && basis(0) == 0 && basis(1) == 1)) {
							if (m2 == -1) {
								std::cout << "Top Generator verified 2*u2sigma^" << (n2) << " / ulambda^" << (-m2) << "\n";
							}
							else {
								std::cout << "Top Generator verified 4*u2sigma^" << (n2) << " / ulambda^" << (-m2) << "\n";

							}
						}
						else {
							throw(8);
						}
					}
					else {
						if ((basis.size() == 1 && basis(0) == 1) || (basis.size() == 2 && (basis(0) == 0 || basis(0) == 1) && basis(1) == 1)) {
							std::cout << "Top Generator verified 4*u2sigma^" << (n2) << " / ulambda^" << (-m2) << "\n";
						}
						else {
							throw(8);
						}
					}
				}
			}
		}


		for (int n1 = 0; n1 <= 1; n1++) {
			for (int n2 = 1; n2 <= rangeN2; n2++) {
				for (int m1 = -1; m1 >= -rangeM1; m1--) {
					for (int m2 = -1; m2 >= -rangeM2; m2--) {
						element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda + s3;
						if (m1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
						}
						else if (m2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
						}
						if (!silence) {

							if (basis.size() == 1 && basis(0) == 1) {
								std::cout << "Top Generator verified asigma^" << n1 << "*u2sigma^" << n2 << "*s3 / alambda^" << (-m1) << "*ulambda^" << (-m2) << "\n";
							}
							else {
								throw(9);
							}
						}
					}
				}
			}
		}

		for (int n1 = 3; n1 <= rangeN1; n1++) {
			for (int n2 = 1; n2 <= rangeN2; n2++) {
				for (int m1 = -1; m1 >= -rangeM1; m1--) {
					element = n1 * asigma + n2 * u2sigma + m1 * alambda;
					basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
					if (!silence) {

						if ((basis.size() == 1 && basis(0) == 1) || (basis.size() == 2 && basis(0) == 1 && basis(1) == 0)) {
							std::cout << "Top Generator verified asigma^" << n1 << "*u2sigma^" << n2 << " / alambda^" << (-m1) << "\n";
						}
						else {
							throw(10);
						}
					}
				}
			}
		}
	}
}
=======
#pragma once
#include <Mackey/Compute.h>
#include <string>
#include <iostream>

///@file
///@brief Contains testing functions for the G=C4 case.
namespace {

	std::string C4MackeyAnswer(int k, int n, int m) {
		std::string MackeyName;
		std::vector<int> sphere = { n,m };
		k = Mackey::invReindex(k, sphere);
		if (n >= 0 && m >= 0) {
			if (k == n + 2 * m && abs(n) % 2 == 0)
				MackeyName = "Z";
			else if (k == n + 2 * m && abs(n) % 2 == 1)
				MackeyName = "Z_-";
			else if (((abs(n) % 2 == 0 && k < n) || (abs(n) % 2 == 1 && k < n + 2 * m)) && abs(k) % 2 == 0 && k >= 0)
				MackeyName = "Z/2";
			else if (abs(n) % 2 == 0 && k >= n && k < n + 2 * m && abs(k) % 2 == 0)
				MackeyName = "Z/4";
			else if (abs(n) % 2 == 1 && k >= n && k < n + 2 * m && abs(k) % 2 == 1)
				MackeyName = "overline Z/2";
			else
				MackeyName = "0";
		}
		else if (n <= 0 && m <= 0) {
			if (k == n + 2 * m && n == -1 && m == 0)
				MackeyName = "Z_-";
			else if (k == n + 2 * m && m == 0 && n != 0 && !(abs(n) % 2))
				MackeyName = "p^*L";
			else if (k == n + 2 * m && n <= -3 && m == 0 && (abs(n) % 2) == 1)
				MackeyName = "p^*L_-";
			else if (k == n + 2 * m && m != 0 && (abs(n) % 2) == 0)
				MackeyName = "L";
			else if (k == n + 2 * m && m != 0 && (abs(n) % 2) == 1)
				MackeyName = "L_-";
			else if (n < 0 && k <= 0 && abs(k) >= 3 && abs(k) % 2 == 1 && ((m == 0 && abs(k) < abs(n)) || (m != 0 && (abs(n) % 2) == 0 && abs(k) <= abs(n) + 1) || (m != 0 && (abs(n) % 2) == 1 && abs(k) < abs(n) + abs(2 * m))))
				MackeyName = "Z/2";
			else if (abs(n) % 2 == 0 && k <= 0 && abs(k) >= abs(n) + 3 && abs(k) < abs(n + 2 * m) && abs(k) % 2 == 1)
				MackeyName = "Z/4";
			else if (abs(n) % 2 == 1 && k <= 0 && abs(k) >= abs(n) + 3 && abs(k) < abs(n + 2 * m) && abs(k) % 2 == 0)
				MackeyName = "overline Z/2";
			else
				MackeyName = "0";
		}
		else if (n > 0 && m < 0) {
			m = -m;
			if (n % 2 == 1) {
				if ((n - 2 * m < k && k <= n - 4 && abs(k) % 2 == 1) || (0 <= k && k < n - 2 * m && abs(k) % 2 == 0))
					MackeyName = "Z/2";
				else if (k == n - 2 * m && m >= 2)
					MackeyName = "L_-";
				else if (k == n - 2 && m == 1)
					MackeyName = "Z_-^flat";
				else if (k == n - 3 && n >= 3 && m >= 2)
					MackeyName = "Q^sharp";
				else if ((n - 2 * m < k && k <= n - 5 && k < 0 && abs(k) % 2 == 0) || (k == -2 && n == 1 && m >= 2))
					MackeyName = "overline Z/2";
				else if (n - 2 * m < k && k <= n - 5 && 0 <= k && abs(k) % 2 == 0)
					MackeyName = "Z/2+overline Z/2";
				else
					MackeyName = "0";
			}
			else {
				if (0 <= k && k <= n - 4 && abs(k) % 2 == 0 && k != n - 2 * m)
					MackeyName = "Z/2";
				else if (n - 2 * m < k && k < n - 3 && abs(k) % 2 == 1)
					MackeyName = "Z/4";
				else if (k == n - 2 * m && n - 2 * m < 0 && m >= 2)
					MackeyName = "L";
				else if (k == n - 2 && m == 1)
					MackeyName = "L^sharp";
				else if (k == n - 3 && m >= 2)
					MackeyName = "Q^sharp";
				else if (k == n - 2 * m && n - 2 * m >= 0 && m >= 2)
					MackeyName = "Z/2+L";
				else
					MackeyName = "0";
			}
		}
		else {
			n = -n;
			if (n % 2 == 0) {
				if (k == -n + 2 * m)
					MackeyName = "Z";
				else if (-n + 1 <= k && k <= -3 && abs(k) % 2 == 1)
					MackeyName = "Z/2";
				else if (-n + 2 <= k && k < -n + 2 * m && abs(k) % 2 == 0)
					MackeyName = "Z/4";
				else if (k == -n)
					MackeyName = "Q";
				else
					MackeyName = "0";
			}
			else {
				if (k == -n + 2 * m && k >= -1)
					MackeyName = "Z_-";
				else if ((-n + 1 <= k && k < 2 * m - n && abs(k) % 2 == 0) || (2 * m - n < k && k <= -3 && abs(k) % 2 == 1 && k != 2 * m - n))
					MackeyName = "Z/2";
				else if (-1 <= k && k < -n + 2 * m && abs(k) % 2 == 1)
					MackeyName = "overline Z/2";

				else if (k == -n && k <= -3)
					MackeyName = "Q";
				else if (k >= -n + 2 && k < -n + 2 * m && k <= -3 && abs(k) % 2 == 1)
					MackeyName = "Z/2+overline Z/2";

				else if (k == -n + 2 * m && k <= -3)
					MackeyName = "Z/2+Z_-";
				else
					MackeyName = "0";
			}
		}
		return MackeyName;
	}

}

///The namespace for testing the answers produced by Mackey (in case of G=C4).
namespace C4Test{

	using namespace Mackey;

	///Tests the additive structure in the range rangeNLow<=n<=rangeNHigh and rangeMLow<=m<=rangeMHigh. silence=1 if we don't want console output.
	template<typename rank_t, typename diff_t>
	void C4MackeyTest(int rangeNLow, int rangeNHigh, int rangeMLow, int rangeMHigh, bool silence)
	{
		std::vector<int> minimum = { rangeNLow, rangeMLow };
		std::vector<int> maximum = { rangeNHigh, rangeMHigh };
		auto degrees = DegreeConstruction(minimum, maximum);
		for (const auto& i : degrees)
		{
			auto n = i[0];
			auto m = i[1];
			auto M = ROHomology<rank_t, diff_t>(i);
			int k = 0;
			for (auto& j : M) {
				if (!silence) {
					auto name = identify(j);
					auto MackeyName = C4MackeyAnswer(k, n, m);
					if (MackeyName != name)
						throw(0);
					else {
					std::cout << "The k=" << invReindex<std::vector<int>>(k, { n,m }) << " homology of the n=" << n << " and m=" << m << " sphere is " << name << "\n";
					}
				}
				k++;
			}
		}
	}



	///Tests the additive structure using powers of asigma, u2sigma, alambda, ulambda in a symmetric range described by the rangeN1,rangeN2,rangeM1,rangeM2 and rangeMLow<=m<=rangeMHigh respectively. silence=1 if we don't want console output.
	template<typename rank_t, typename diff_t>
	void C4Multtest(int rangeN1, int rangeN2, int rangeM1, int rangeM2, bool silence) {
		rank_t basis, element, elementnow, elementless;
		rank_t asigma(3), u2sigma(3), alambda(3), ulambda(3), w3(3), x11(3), s3(3);
		asigma << 0, 1, 0; u2sigma << 2, 2, 0; alambda << 0, 0, 1; ulambda << 2, 0, 1; w3 << -3, -3, 0; x11 << -3, -1, -1; s3 << -3, 0, -2;


		for (int n1 = 0; n1 <= rangeN1; n1++) {
			for (int n2 = 0; n2 <= rangeN2; n2++) {
				for (int m1 = 0; m1 <= rangeM1; m1++) {
					for (int m2 = 0; m2 <= rangeM2; m2++) {
						if (n1 <= 1 || m2 == 0) {
							elementnow = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda;
							if (n1 == 0 && n2 == 0 && m1 == 0 && m2 == 0) {
								continue;
							}
							if (n1 > 0) {
								elementless = elementnow - asigma;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, elementless);
							}
							else if (n2 > 0) {
								elementless = elementnow - u2sigma;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, elementless);
							}
							else if (m1 > 0) {
								elementless = elementnow - alambda;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, elementless);
							}
							else if (m2 > 0) {
								elementless = elementnow - ulambda;
								basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, elementless);
							}
							if (!silence) {
								if (basis.size() == 1 && basis(0) == 1) {
									std::cout << "Top Generator verified asigma^" << n1 << "*u2sigma^" << n2 << "*alambda^" << m1 << "*ulambda^" << m2 << "\n";
								}
								else {
									throw(1);
								}
							}
						}
					}
				}
			}
		}


		basis = ROGreen<rank_t, diff_t, rank_t>(2, 2 * asigma, ulambda);
		if (!silence) {
			if (basis.size() == 1 && basis(0) == 2) {
				std::cout << "Gold Verified!" << "\n";

			}
			else {
				throw(2);
			}
		}


		for (int n2 = 0; n2 >= -rangeN2; n2--) {
			for (int m2 = 0; m2 >= -rangeM2; m2--) {
				element = n2 * u2sigma + m2 * ulambda;
				if (n2 == 0 && m2 == 0) {
					continue;
				}
				else if (n2 < 0) {
					basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
				}
				else {
					basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
				}
				if (!silence) {

					if (n2 == -1 && m2 == 0) {
						if (basis.size() == 1 && basis(0) == 2) {
							std::cout << "Top Generator verified: 2 / u2sigma^" << (-n2) << "\n";
						}
						else {
							throw(3);
						}
					}
					else if (n2 == 0 && m2 == -1) {
						if (basis.size() == 1 && basis(0) == 4) {
							std::cout << "Top Generator verified: 4 / ulambda^" << (-m2) << "\n";
						}
						else {
							throw(3);
						}
					}
					else {
						if (basis.size() == 1 && basis(0) == 1) {
							if (m2 == 0) {
								std::cout << "Top Generator verified: 2 / u2sigma^" << (-n2) << "\n";
							}
							else {
								std::cout << "Top Generator verified: 4 / u2sigma^" << (-n2) << "*ulambda^" << (-m2) << "\n";
							}
						}
						else {
							throw(3);
						}
					}
				}
			}
		}

		for (int n2 = 0; n2 >= -rangeN2; n2--) {
			for (int m1 = 0; m1 >= -rangeM1; m1--) {
				for (int m2 = 0; m2 >= -rangeM2; m2--) {
					element = n2 * u2sigma + m1 * alambda + m2 * ulambda + s3;
					if (n2 == 0 && m1 == 0 && m2 == 0) {
						continue;
					}
					if (n2 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
					}
					else if (m1 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
					}
					else if (m2 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
					}
					if (!silence) {

						if (basis.size() == 1 && basis(0) == 1) {
							std::cout << "Top Generator verified s3 / u2sigma^" << (-n2) << "*alambda^" << (-m1) << "*ulambda^" << (-m2) << "\n";
						}
						else {
							throw(4);
						}
					}
				}
			}
		}

		for (int n1 = 0; n1 >= -rangeN1; n1--) {
			for (int n2 = 0; n2 >= -rangeN2; n2--) {
				for (int m1 = 0; m1 >= -rangeM1; m1--) {
					for (int m2 = 0; m2 >= -rangeM2; m2--) {
						if (n1 < 0 && m2 < 0) {
							continue;
						}
						if (n1 == 0 && n2 == 0 && m1 == 0 && m2 == 0) {
							continue;
						}
						element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda + x11;
						if (n1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, element);
						}
						if (n2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
						}
						else if (m1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
						}
						else if (m2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
						}
						if (!silence) {

							if (basis.size() == 1 && basis(0) == 1) {
								std::cout << "Top Generator verified x11 / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "*alambda^" << (-m1) << "*ulambda^" << (-m2) << "\n";
							}
							else {
								throw(5);
							}
						}
					}
				}
			}
		}

		for (int n2 = -1; n2 >= -rangeN2; n2--) {
			for (int m2 = 1; m2 <= rangeM2; m2++) {
				element = n2 * u2sigma + m2 * ulambda;
				basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
				if (!silence) {

					if (basis.size() == 1 && basis(0) == 1) {
						std::cout << "Top Generator verified ulambda^" << (m2) << " / u2sigma^" << (-n2) << "\n";
					}
					else {
						throw(1);
					}
				}
			}
		}



		for (int n1 = 0; n1 >= -1; n1--) {
			for (int n2 = 0; n2 >= -rangeN2; n2--) {
				for (int m1 = 1; m1 <= rangeM1; m1++) {
					for (int m2 = 0; m2 <= rangeM2; m2++) {
						element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda;
						if (n1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, element);
						}
						else if (n2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
						}
						if (!silence) {

							if ((n1 == -1 && n2 == 0 && m2 == 0) || (n1 == -1 && m2 != 0) || (n1 == 0 && n2 == -1 && m2 == 0)) {
								if (basis.size() == 1 && basis(0) == 2) {
									if (m2 == 0 || n1 != 0) {
										std::cout << "Top Generator verified 2*alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";
									}
									else {
										std::cout << "Top Generator verified alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";
									}
								}
								else {
									throw(6);
								}
							}
							else {
								if (basis.size() == 1 && basis(0) == 1) {
									if (m2 == 0 || n1 != 0) {
										std::cout << "Top Generator verified 2*alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << " \n";
									}
									else {
										std::cout << "Top Generator verified alambda^" << m1 << "*ulambda^" << m2 << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";

									}
								}
								else {
									throw(6);
								}
							}
						}
					}
				}
			}
		}



		for (int n1 = 0; n1 >= -rangeN1; n1--) {
			for (int n2 = 0; n2 >= -rangeN2; n2--) {
				for (int m1 = 1; m1 <= rangeM1; m1++) {
					element = n1 * asigma + n2 * u2sigma + m1 * alambda + w3;
					if (n1 == 0 && n2 == 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, m1 * alambda, w3);
					}
					if (n1 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, asigma, element);
					}
					else if (n2 < 0) {
						basis = ROGreen<rank_t, diff_t, rank_t>(2, u2sigma, element);
					}
					if (!silence) {

						if (basis.size() == 1 && basis(0) == 1) {
							std::cout << "Top Generator verified alambda^" << m1 << "*w3" << " / asigma^" << (-n1) << "*u2sigma^" << (-n2) << "\n";

						}
						else {
							throw(7);
						}
					}
				}
			}
		}

		for (int n2 = 1; n2 <= rangeN2; n2++) {
			for (int m2 = -1; m2 >= -rangeM2; m2--) {
				element = n2 * u2sigma + m2 * ulambda;
				basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element, 0, 1);
				if (!silence) {

					if (m2 == -2 || m2 == -1) {
						if ((basis.size() == 1 && basis(0) == 2) || (basis.size() == 2 && basis(0) == 0 && basis(1) == 1)) {
							if (m2 == -1) {
								std::cout << "Top Generator verified 2*u2sigma^" << (n2) << " / ulambda^" << (-m2) << "\n";
							}
							else {
								std::cout << "Top Generator verified 4*u2sigma^" << (n2) << " / ulambda^" << (-m2) << "\n";

							}
						}
						else {
							throw(8);
						}
					}
					else {
						if ((basis.size() == 1 && basis(0) == 1) || (basis.size() == 2 && (basis(0) == 0 || basis(0) == 1) && basis(1) == 1)) {
							std::cout << "Top Generator verified 4*u2sigma^" << (n2) << " / ulambda^" << (-m2) << "\n";
						}
						else {
							throw(8);
						}
					}
				}
			}
		}


		for (int n1 = 0; n1 <= 1; n1++) {
			for (int n2 = 1; n2 <= rangeN2; n2++) {
				for (int m1 = -1; m1 >= -rangeM1; m1--) {
					for (int m2 = -1; m2 >= -rangeM2; m2--) {
						element = n1 * asigma + n2 * u2sigma + m1 * alambda + m2 * ulambda + s3;
						if (m1 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
						}
						else if (m2 < 0) {
							basis = ROGreen<rank_t, diff_t, rank_t>(2, ulambda, element);
						}
						if (!silence) {

							if (basis.size() == 1 && basis(0) == 1) {
								std::cout << "Top Generator verified asigma^" << n1 << "*u2sigma^" << n2 << "*s3 / alambda^" << (-m1) << "*ulambda^" << (-m2) << "\n";
							}
							else {
								throw(9);
							}
						}
					}
				}
			}
		}

		for (int n1 = 3; n1 <= rangeN1; n1++) {
			for (int n2 = 1; n2 <= rangeN2; n2++) {
				for (int m1 = -1; m1 >= -rangeM1; m1--) {
					element = n1 * asigma + n2 * u2sigma + m1 * alambda;
					basis = ROGreen<rank_t, diff_t, rank_t>(2, alambda, element);
					if (!silence) {

						if ((basis.size() == 1 && basis(0) == 1) || (basis.size() == 2 && basis(0) == 1 && basis(1) == 0)) {
							std::cout << "Top Generator verified asigma^" << n1 << "*u2sigma^" << n2 << " / alambda^" << (-m1) << "\n";
						}
						else {
							throw(10);
						}
					}
				}
			}
		}
	}
}
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
