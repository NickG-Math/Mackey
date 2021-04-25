#pragma once
#include <iostream>
#include <array>
#include <Eigen/Dense>
#include "Types/SFINAE.hpp"

///@file
///@brief Contains the class \ref mackey::Z_mod

namespace mackey {
	///////////////////////////////////
	///The class of Z/N coefficients where N is prime.

	///The operators are self explanatory. User must ensure that N is prime for division to work.
	/////////////////////////////////
	template<int64_t N, typename T = int64_t>
	class Z_mod {
	public:
		constexpr static int64_t order = N; ///<The N of Z/N
		T x; ///< A modulo N number
		Z_mod() : x(0) {} ///<Default value 0
		Z_mod(bool x) : x(x) {} ///<Initialize from 0,1
		Z_mod(int x); ///<Initialize from int
		Z_mod(int64_t x); ///<Initialize from 64bit int
		explicit operator char() const;
		explicit operator short() const;
		explicit operator int() const;
		explicit operator int64_t() const;
		explicit operator unsigned char() const;
		explicit operator unsigned short() const;
		explicit operator unsigned int() const;
		explicit operator uint64_t() const;
		Z_mod<N, T> operator +(Z_mod<N, T> b) const;
		Z_mod<N, T> operator -(Z_mod<N, T> b) const;
		Z_mod<N, T>& operator +=(Z_mod<N, T> b);
		Z_mod<N, T>& operator -=(Z_mod<N, T> b);
		Z_mod<N, T>& operator *=(Z_mod<N, T> b);
		Z_mod<N, T>& operator /=(Z_mod<N, T> b);
		bool operator ==(Z_mod<N, T> a) const;
		bool operator !=(Z_mod<N, T> a) const;
		bool operator <=(Z_mod<N, T> a) const; ///<Needed for Eigen pruning; standard order on 0,...,N-1
	};


	template<int64_t N, typename T>
	Z_mod<N, T> operator -(Z_mod<N, T> a);

	template<int64_t N, typename T>
	Z_mod<N, T> operator *(Z_mod<N, T> a, Z_mod<N, T> b); //Eigen needs this to be non member

	template<int64_t N, typename T>
	Z_mod<N, T> operator /(Z_mod<N, T> a, Z_mod<N, T> b);

	///The usual absolute value for integer and Z/N types
	template<typename T>
	T abs(T a);

	template<int64_t N, typename T> //Eigen needs this to be non member
	std::ostream& operator<<(std::ostream& out, const Z_mod<N, T> a);

	///The \f$\mathbf Z/2\f$ coefficients
	using Z2 = Z_mod<2, bool>;
}

///See the Eigen documentation for this. 
namespace Eigen {
	using namespace mackey;

	///Specializing NumTraits to Z/nZ coefficients
	template<int64_t N, typename T>
	struct NumTraits<Z_mod<N, T>>
	{
		typedef Z_mod<N, T> Real;
		typedef Z_mod<N, T> Nested;
		typedef int Literal;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 0,
			RequireInitialization = 0,
			ReadCost = 1,
			AddCost = 1,
			MulCost = 1
		};
		static inline Z_mod<N, T> dummy_precision() { return Z_mod<N, T>(0); }
		static inline int digits10() { return 0; }
	};
}

#include "impl/Z_n.ipp"
