#pragma once
#include <vector>
#include "General.h"

///@file
///@brief Contains the class for formal polynomials in arbitrary coefficients


namespace Mackey {
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
		return Polynomial<T>(-a.p);
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
}
