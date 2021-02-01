#include "Wrapper.h"


//Demo for a group specific implementation for C_2^n. Currently for n=2, you only have to set the power variable below to the n you want

//First set the variables

const int GroupSpecific::Variables::prime = 2;
const int GroupSpecific::Variables::power = 2;
const int GroupSpecific::Variables::reps = GroupSpecific::Variables::power;

std::vector<int> power2_dimensions(int power) {
	std::vector<int> dims(power, 2);
	dims[0] = 1;
	return dims;
}

const std::vector<int> GroupSpecific::Variables::sphere_dimensions = power2_dimensions(GroupSpecific::Variables::power);

//Anonymous namespace to setup the recursion used in obtaining the positive chains for any C_2^n. 
namespace {
	///////////////////////////////////
	///Formal polynomials with arbitrary coefficients

	///Currently only used for a recursion in C2n_Implementation. The operators are self explanatory.
	/////////////////////////////////
	template<typename T>
	class Polynomial {
	public:
		std::vector<T> p; ///<Polynomial v_0+...+v_nx^n stored as vector (v_0,...,v_n)
		int degree; ///<Degree of the polynomial

		///Default value 0
		Polynomial() : p({ 0 }) { degree = -1; }

		///Construct from vector v_0,...,v_n the polynomial v_0+...+v_nx^n
		Polynomial(const std::vector<T>& q) : p(q) {
			while (!p.empty() && p.back() == 0) {
				p.pop_back();
			}
			degree = (int)p.size() - 1; //avoid usigned underflow
		}

		///Construct the monomial ax^n given a,n
		Polynomial(T coeff, int power) {
			if (coeff != 0) {
				p.resize(power + 1); p[power] = coeff; degree = power;
			}
			else {
				p = { 0 };
				degree = -1;
			}
		}

		///Construct the constant a
		Polynomial(T coeff) :Polynomial<T>(coeff, 0) {}

		///Pad polynomial with extra zeros after the highest power
		void pad(int m) {
			p.resize(m);
		}
	};

	template<typename T>
	Polynomial<T> operator +(const Polynomial<T>& a, const Polynomial<T>& b) {
		std::vector<T> sum;
		sum.reserve(std::max(a.degree, b.degree) + 1);
		for (int i = 0; i < std::max(a.degree, b.degree) + 1; i++) {
			T summand = 0;
			if (i <= a.degree)
				summand += a.p[i];
			if (i <= b.degree)
				summand += b.p[i];
			sum.push_back(summand);
		}
		return Polynomial(sum);
	}

	template<typename T>
	Polynomial<T>& operator +=(Polynomial<T>& a, const Polynomial<T>& b) {
		a = a + b;
		return a;
	}


	template<typename T>
	Polynomial<T> operator -(const Polynomial<T>& a) {
		return Polynomial<T>(Mackey::operator-(a.p));
	}

	template<typename T>
	Polynomial<T> operator -(const Polynomial<T>& a, const Polynomial<T>& b) {
		return (a + (-b));
	}

	template<typename T>
	Polynomial<T>& operator -=(Polynomial<T>& a, const Polynomial<T>& b) {
		a = a - b;
		return a;
	}

	template<typename T>
	Polynomial<T> operator *(const Polynomial<T>& a, const Polynomial<T>& b) {
		std::vector<T> product;
		product.reserve(a.degree + b.degree + 1);
		for (int i = 0; i < a.degree + b.degree + 1; i++) {
			T sum = 0;
			int limit_low = std::max(0, (int)(i - b.degree));
			int limit_high = std::min(i, (int)(a.degree));
			for (int j = limit_low; j <= limit_high; j++) {
				sum += a.p[j] * b.p[i - j];
			}
			product.push_back(sum);
		}
		return Polynomial(product);
	}


	template<typename T>
	std::pair<Polynomial<T>, Polynomial<T>> divide(const Polynomial<T>& a, const Polynomial<T>& b) {
		Polynomial<T> quotient;
		if (a.degree < b.degree) {
			return std::make_pair(quotient, a);
		}
		Polynomial<T> remainder = a;
		do {
			Polynomial<T> poly(remainder.p.back() / b.p.back(), remainder.degree - b.degree);
			quotient += poly;
			remainder -= poly * b;
		} while (remainder.degree >= b.degree && remainder.degree != -1);

		return std::make_pair(quotient, remainder);
	}


	template<typename T>
	Polynomial<T> operator /(const Polynomial<T>& a, const Polynomial<T>& b) {
		return divide(a, b).first;
	}

	template<typename T>
	Polynomial<T> operator %(const Polynomial<T>& a, const Polynomial<T>& b) {
		return divide(a, b).second;
	}

	///Geometric sum 1+...+x^{max}
	template<typename T>
	Polynomial<T> geometric(int max) {
		std::vector<T> p(max + 1, 1);
		return Polynomial<T>(p);
	}

	//The d and x sequences are defined in HHR17 pg 392.

	template<typename deg_t>
	std::vector<int> dsequence(int n, const deg_t& sphere) {
		std::vector<int> d(sphere.size() + 1);
		for (int i = 0; i < sphere.size(); i++) {
			d[i + 1] = d[i] + sphere[i] * GroupSpecific::Variables::sphere_dimensions[i];
		}
		return d;
	}

	template< typename Scalar, typename deg_t>
	std::vector<Polynomial<Scalar>> xsequence(int n, const deg_t& sphere, const std::vector<int>& d) {
		std::vector<Polynomial<Scalar>> x(n + 1);
		if (n < 2)
			return x;
		x[2] = Polynomial<Scalar>(1);
		for (int i = 3; i <= n; i++) {
			int j = 1;
			while (2 + d[j - 1] >= i || i > 2 + d[j])
				j++;
			x[i] = geometric<Scalar>((1 << j) - 1) / x[i - 1];
		}
		return x;
	}

	using namespace Mackey;

	//Then write the standard chains. The function signature and name need not be changed, only the definition
	template<typename rank_t, typename diff_t, typename deg_t>
	Chains<rank_t, diff_t> myfunction(int n, const deg_t& sphere) {

		auto d = dsequence(n, sphere);
		auto x = xsequence<typename diff_t::Scalar>(n, sphere, d);
		std::vector<rank_t> rank;
		std::vector<diff_t> diff;

		rank.reserve(n + 1);
		rank.push_back(altmatrix<rank_t>(1, 1, { 1 }));
		for (int i = 1; i <= n; i++) {
			int j = 1;
			while (d[j - 1] >= i || i > d[j])
				j++;
			rank.push_back(altmatrix<rank_t>(1, 1, { static_cast<typename rank_t::Scalar>(1 << j) }));
		}
		diff.resize(n + 1);
		if (rank.size() == 1) {
			return Chains<rank_t, diff_t>(rank, diff);
		}
		else {
			diff[1] = altmatrix<diff_t>(rank[0][0], rank[1][0], { 1 });
			std::vector<typename diff_t::Scalar> v = { 1,-1 };
			Polynomial<typename diff_t::Scalar> r(v); //1-g polynomial
			for (int i = 2; i <= n; i++) {
				auto relation = Polynomial<typename diff_t::Scalar>(1) - Polynomial<typename diff_t::Scalar>(1, rank[i - 1][0]);
				if (i % 2) {//odd
					auto s = x[i] % relation;
					s.pad(rank[i - 1][0]);
					diff[i] = altmatrix<diff_t>(rank[i - 1][0], rank[i][0], s.p);
				}
				else { //even
					auto rxd = (r * x[i]) % relation;
					rxd.pad(rank[i - 1][0]);
					diff[i] = altmatrix<diff_t>(rank[i - 1][0], rank[i][0], rxd.p);
				}
			}
			return Chains<rank_t, diff_t>(rank, diff);
		}
	}
}


//Finally set the standard chains. No need to change the syntax in what follows:
template<typename rank_t, typename diff_t, typename deg_t>
const typename GroupSpecific::Function<rank_t, diff_t, deg_t>::functype GroupSpecific::Function<rank_t, diff_t, deg_t>::PositiveChains = &myfunction;
