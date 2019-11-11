//Eigen directives
#define EIGEN_USE_MKL_ALL
//#define EIGEN_NO_DEBUG //disable Eigen assertions
//#define NDEBUG

//Mackey directives

#ifdef _MSC_VER
#define _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING //Eigen won't work otherwise.
#endif
//note: These preprocessor directives must be before any includes.

#include <iostream>
#include "C4Verify.h"


using namespace C4Test;

typedef Eigen::Matrix<short, 1, -1> rank_t;
typedef Eigen::Matrix<char, -1, -1> diff_t;



inline std::array<int, 4> interfaceAdd() {
	std::array<int, 4> ranges;
	bool flag = 1;
	while (flag) {
		for (int i = 0; i < 4; i++) {
			std::cin >> ranges[i];
		}
		if (ranges[0] > ranges[1] || ranges[2] > ranges[3]) {
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

inline std::array<int, 4> interfaceMul() {
	std::array<int, 4> ranges;
	bool flag = 1;
	while (flag) {
		for (int i = 0; i < 4; i++) {
			std::cin >> ranges[i];
		}
		if (ranges[0] < 0 || ranges[1] < 0 || ranges[2] < 0 || ranges[3] < 0) {
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



inline void showandtell() {
	std::cout << "The RO(C4) homology of a point in Z coefficients \n\n";
	std::cout << "Please enter \"a\" (no quotes) if you want the additive structure and \"m\" (no quotes) if you want the multiplicative one\n";
	bool flag = 1;
	std::string choice;
	while (flag) {
		std::cin >> choice;
		if (choice != "a" && choice != "m") {
			std::cout << "Please enter \"a\" (no quotes) if you want the additive structure and \"m\" (no quotes) if you want the multiplicative one\n";
		}
		else {
			flag = 0;
		}

	}

	if (choice == "a") {
		std::cout << "Additive structure selected.\nTo compute H_kS^{n*sigma+ m*lambda} for a<=n<=b and c<=m<=d enter the a b c d separated by spaces\nExample: -3 0 -5 10\n";
		auto ranges = interfaceAdd();
		std::cin.get();
		C4MackeyTest<rank_t, diff_t>(ranges[0], ranges[1], ranges[2], ranges[3], 0);
	}
	if (choice == "m") {
		std::cout << "Multiplicative structure selected.\nWe will write all generators in H_kS^{n*sigma+ m*lambda} in terms of Euler classes asigma, u2sigma, orientation classes alambda, ulambda, the transfers w3 and x11 and finally the generator s3\n";
		std::cout << "Please enter the maximum powers the asigma, u2sigma, alambda, ulambda can take (either in numerators or denominators) as a b c d respectively, separated by spaces.\nExample 3 2 1 5\n";
		auto ranges = interfaceMul();
		std::cin.get();
		C4Multtest<rank_t, diff_t>(ranges[0], ranges[1], ranges[2], ranges[3], 0);
	}

	std::cout << "\n\n Computation successful. Press any key to exit.";
	std::cin.get();
}



int main() {
	showandtell();
	return 0;
}
