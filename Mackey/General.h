#pragma once
#include <vector>
#include <Eigen/Dense>
#include <numeric> //lcm

///@file
///@brief Contains general functions and vector overloads independent of everything else.

///The namespace for everything that is group-independent.
namespace Mackey {

	///////////////////////////////////////////////////
	///Given minimum and maximum degrees, constructs everything in between.

	///deg_t here is usually an std::vector<int>.
	///For example if minimum=[-1,-1] and maximum=[2,3] then the result is {[-1,-1],[0,-1],[1,-1],[2,-1],[-1,0],...,[2,3]}
	///////////////////////////////////////////////////
	template<typename deg_t>
	std::vector<deg_t> DegreeConstruction(const deg_t& minimum, const deg_t& maximum) {
		auto totallength = maximum[0] - minimum[0] + 1;
		for (decltype(maximum.size()) i = 1; i < maximum.size(); i++) {
			totallength *= maximum[i] - minimum[i] + 1;
		}
		std::vector<deg_t> degrees;
		auto degree = minimum;
		degrees.reserve(totallength);
		for (int counter = 0; counter < totallength; counter++) {
			degrees.push_back(degree);
			decltype(maximum.size()) i = 0;
			while (i < maximum.size() - 1 && degree[i] == maximum[i]) {
				degree[i] = minimum[i];
				i++;
			}
			degree[i]++;
		}
		return degrees;
	}


	/////////////////////////////////////////////
///Given vector or array returns vector starting from index start.

///Example given v apply as tail(v.data(),v.size(), start)
////////////////////////////////////////////
	template<typename T>
	inline std::vector<T> tail(const T* const& ptr, int size, int start) {
		std::vector<T> vec(ptr + start, ptr + size);
		return vec;
	}


	/////////////////////////////////////////////
	///Given vector or array returns vector starting from index 1.

	///Example given v apply as tail(v.data(),v.size())
	////////////////////////////////////////////
	template<typename T>
	inline std::vector<T> tail(const T* const& ptr, int size) {
		return tail(ptr,size,1);
	}

	///Find first instance of a in v
	template<typename T>
	int find(const T& v, int a) {
		for (decltype(v.size()) i = 0; i < v.size(); i++) {
			if (v[i]==a)
				return i;
		}
		return -1;
	}

	///Count the instances of a in v and return the last one
	template<typename T, typename S>
	std::pair<int,int> findcount(const T& v, S a) {
		int counter = 0;
		int pos = 0;
		for (int i = 0; i < v.size(); i++) {
			if (v[i] == a) {
				counter++;
				pos = i;
			}
		}
		return std::make_pair(counter, pos);
	}

	///Makes the multiple of the basis vector 0,...,multiple,...,0
	inline std::vector<int> basisVector(int total, int position, int multiple) {
		std::vector<int> v(total);
		v[position] = multiple;
		return v;
	}

	///Makes the basis vector 0,...,1,...,0
	inline std::vector<int> basisVector(int total, int position) {
		return basisVector(total,position,1);
	}



	///Makes a new matrix out of A keeping only its rows indicated by Z
	template<typename Derived, typename T>
	Derived KeepRow(const Eigen::DenseBase<Derived>& A, const std::vector<T>& Z) { //when the new stable version of Eigen releases this will be deprecated
		Derived B(Z.size(), A.cols());
		int j = 0;
		for (auto const& i : Z) {
			B.row(j) = A.row(i);
			j++;
		}
		return B;
	}

	///Makes a new matrix out of A keeping only its columns indicated by Z
	template<typename Derived, typename T>
	Derived KeepCol(const Eigen::DenseBase<Derived>& A, const std::vector<T>& Z) { //when the new stable version of Eigen releases this will be deprecated
		Derived B(A.rows(), Z.size());
		int j = 0;
		for (auto const& i : Z) {
			B.col(j) = A.col(i);
			j++;
		}
		return B;
	}

	////////////////////////////////////////////////////////////
	///Takes the sum of the entries of an Eigen vector at higher precision than its inputs.
	template<typename Derived>
	int summation(const Eigen::MatrixBase<Derived>& A) {
		return A.template cast<int>().sum();
	}

	///Circularly rotates given V (an Eigen/std vector)
	template<typename T>
	void rotate(T& V) {
		for (int i = 0; i < V.size() - 1; i++) {
			std::swap(V[i], V[i + 1]);
		}
	}

	////////////////////////////////////////////////////////////////
	///Creates a matrix of size m x n that alternates the given vector v

	///Eg for m=2, n=2, v={a,b} we get {a,b // b,a}
	//////////////////////////////////////////////////////////////
	template<typename Matrix_t>
	Matrix_t altmatrix(int m, int n, const std::vector<typename Matrix_t::Scalar>& v) {
		typedef Eigen::Matrix<typename Matrix_t::Scalar, -1, 1> vector_t;
		Matrix_t matrix(m, n);
		vector_t alt = Eigen::Map<const vector_t>(v.data(), v.size());
		vector_t column(m);
		for (int i = 0; i < m; i += alt.size()) {
			column.segment(i, alt.size()) = alt;
		}
		for (int j = n-1; j>=0; j--) {
			rotate(column);
			matrix.col(j)=column;
		}
		return matrix;
	}
	   	  

	///Coordinate-wise sum of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(std::min(a.size(), b.size()));
		for (size_t i = 0; i < std::min(a.size(), b.size()); i++) {
			c.push_back(a[i] + b[i]);
		}
		return c;
	}

	///Coordinate-wise difference of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(std::max(a.size(), b.size()));
		for (size_t i = 0; i < std::min(a.size(), b.size()); i++) {
			c.push_back(a[i] - b[i]);
		}
		return c;
	}

	///Coordinate-wise opposite of vector
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a)
	{
		std::vector<T> b;
		b.reserve(a.size());
		for (const auto & i:a) {
			b.push_back(-i);
		}
		return b;
	}

	///Coordinate-wise difference of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator*(T a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(b.size());
		for (int i = 0; i < b.size(); i++) {
			c.push_back(a * b[i]);
		}
		return c;
	}

	///Least common multiple of vector of elements
	template<typename T>
	int lcm(const T& v) {
		if (v.size()==0)
			return 0;
		auto n=1;
		for (const auto& i:v) {
			n = std::lcm(n, i);
		}
		return n;
	}



}
