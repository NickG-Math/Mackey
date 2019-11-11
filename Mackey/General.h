#pragma once
#include <vector>
#include <Eigen/Dense>

///@file
///@brief Contains general functions and vector overloads independent of everything else.

///The namespace for everything that is group-independent.
namespace Mackey {

	///////////////////////////////////////////////////
	///Given minimum and maximum degrees, constructs everything in between.
	///
	///deg_t here is usually an std::vector<int>.
	///For example if minimum=[-1,-1] and maximum=[2,3] then the result is {[-1,-1],[0,-1],[1,-1],[2,-1],[-1,0],...,[2,3]}
	///////////////////////////////////////////////////
	template<typename deg_t>
	std::vector<deg_t> DegreeConstruction(const deg_t& minimum, const deg_t& maximum) {
		auto totallength = maximum[0] - minimum[0] + 1;
		for (int i = 1; i < maximum.size(); i++) {
			totallength *= maximum[i] - minimum[i] + 1;
		}
		std::vector<deg_t> degrees;
		auto degree = minimum;
		degrees.reserve(totallength);
		int counter = 0;
		for (int counter = 0; counter < totallength; counter++) {
			degrees.push_back(degree);
			int i = 0;
			while (i < maximum.size() - 1 && degree[i] == maximum[i]) {
				degree[i] = minimum[i];
				i++;
			}
			degree[i]++;
		}
		return degrees;
	}

	/////////////////////////////////////////////
	///Given vector or array returns vector starting from index 1.
	
	///Example given v apply as remove0th(v.data(),v.size())
	////////////////////////////////////////////
	template<typename T>
	inline std::vector<T> remove0th(const T* const& ptr, int size) {
		std::vector<T> vec(ptr + 1, ptr + size);
		return vec;
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

	///T is anything with sum() and cast<int>() methods.
	///////////////////////////////////////////////////////////
	template<typename T>
	int summation(const T& A) {
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

}
