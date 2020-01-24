#pragma once
#include <vector>
#include "SFINAE.h"
#include <numeric> //lcm

///@file
///@brief Contains general functions and vector overloads independent of everything else.

///The namespace for everything that is group-independent.
namespace Mackey {

	///////////////////////////////////////////////////
	///Given minimum and maximum degrees, constructs everything in between.

	///For example if minimum={-1,-1} and maximum={2,3} then the result is {{-1,-1},{0,-1},{1,-1},{2,-1},{-1,0},...,{2,3}}
	///////////////////////////////////////////////////
	template<typename deg_t = std::vector<int>>
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

	///Given a vector of vectors v form all combinations of the form {v[0][?],v[1][?],...} for varying ?
	template<typename T>
	std::vector<std::vector<T>> combinations(const std::vector<std::vector<T>>& v) {
		std::vector<std::vector<T>> w;
		std::vector<T> each;
		each.reserve(v.size());
		std::vector<int> counter(v.size(), 0);
		int totalsize = 1;
		for (const auto& i : v)
			totalsize *= i.size();
		w.reserve(totalsize);
		int change = v.size() - 1;
		while (change >= 0) {
			while (counter[change] == v[change].size()) {
				counter[change] = 0;
				change--;
				if (change < 0)
					return w;
				counter[change]++;
			}
			each.clear();
			for (int i = 0; i < v.size(); i++)
				each.push_back(v[i][counter[i]]);
			w.push_back(each);
			change = v.size() - 1;
			counter[change]++;
		}
		return w;
	}

	/////////////////////////////////////////////
	///Given vector or array returns vector starting from index start.

	///Example: Given v apply as tail(v.data(),v.size(), start)
	////////////////////////////////////////////
	template<typename T>
	inline std::vector<T> tail(const T* const& ptr, int size, int start) {
		return std::vector<T>(ptr + start, ptr + size);
	}

	///Find first instance of a in v for non Eigen matrices
	template<typename T, typename S, std::enable_if_t<!SFINAE::is_Dense<S>::value, int> = 0>
	int find(const T& v, const S& a) {
		for (int i = 0; i < v.size(); i++) {
			if (v[i] == a)
				return i;
		}
		return -1;
	}


	///Find first instance of a in v for dense Eigen matrices
	template<typename T, typename S, std::enable_if_t<SFINAE::is_Dense<S>::value, int> = 0>
	int find(const T& v, const S& a) {
		for (int i = 0; i < v.size(); i++) {
			if (v[i].rows() == a.rows() && v[i].cols() == a.cols() && v[i] == a) //otherwise asserton fails
				return i;
		}
		return -1;
	}

	///Makes the multiple of the basis vector 0,...,multiple,...,0
	template<typename rank_t>
	inline rank_t basisElement(int length, int position, int multiple) {
		rank_t a = rank_t::Zero(length);
		a[position] = multiple;
		return a;
	}

	///Makes the the basis vector 0,...,1,...,0
	template<typename rank_t>
	inline rank_t basisElement(int length, int position) {
		return basisElement<rank_t>(length, position, 1);
	}

	///Makes a new matrix out of A keeping only its rows indicated by Z
	template<typename T, typename S>
	T KeepRow(const T& A, const S& Z) { //when the new stable version of Eigen releases this will be deprecated
		T B(Z.size(), A.cols());
		size_t j = 0;
		for (const auto& i : Z) {
			B.row(j) = A.row(i);
			j++;
		}
		return B;
	}

	///Makes a new matrix out of A keeping only its columns indicated by Z
	template<typename T, typename S>
	T KeepCol(const T& A, const S& Z) { //when the new stable version of Eigen releases this will be deprecated
		T B(A.rows(), Z.size());
		size_t j = 0;
		for (const auto& i : Z) {
			B.col(j) = A.col(i);
			j++;
		}
		return B;
	}

	////////////////////////////////////////////////////////////
	///Takes the sum of the entries of an Eigen vector at higher precision than its inputs.
	template<typename Derived>
	long summation(const Eigen::MatrixBase<Derived>& A) {
		return A.template cast<long>().sum();
	}


	template<typename rank_t>
	long summation(const rank_t& u, long limit) {
		long sum = 0;
		for (long i = 0; i < limit; i++)
			sum += u[i];
		return sum;
	}


	///Circularly rotates given V (an Eigen/std vector)
	template<typename T>
	void rotate(T& V) {
		for (int i = 0; i < V.size() - 1; i++)
			std::swap(V[i], V[i + 1]);
	}

	////////////////////////////////////////////////////////////////
	///Creates a matrix of size m x n that alternates the given vector v

	///Eg for m=2, n=2, v={a,b} we get {a,b // b,a}
	//////////////////////////////////////////////////////////////
	template<typename Matrix_t, std::enable_if_t<SFINAE::is_Dense<Matrix_t>::value, int> = 0 >
	Matrix_t altmatrix(int m, int n, const std::vector<Scalar_t<Matrix_t>>& v) {
		Matrix_t matrix(m, n);
		col_t<Matrix_t> alt = Eigen::Map<const col_t<Matrix_t>>(v.data(), v.size());
		col_t<Matrix_t> column(m);
		for (int i = 0; i < m; i += alt.size())
			column.segment(i, alt.size()) = alt;
		for (int j = n - 1; j >= 0; j--) {
			rotate(column);
			matrix.col(j) = column;
		}
		return matrix;
	}


	///Sparse version of altmatrix
	template<typename Matrix_t, std::enable_if_t<SFINAE::is_Sparse<Matrix_t>::value, int> = 0 >
	Matrix_t altmatrix(int m, int n, const std::vector<Scalar_t<Matrix_t>>& v) {
		auto dense = altmatrix<mat_t<Matrix_t>>(m, n, v);
		return dense.sparseView();
	}

	///Tests if vector is zero
	template <typename T>
	bool isZero(const std::vector<T>& a)
	{
		for (const auto& i : a) {
			if (i != 0)
				return 0;
		}
		return 1;
	}

	///Coordinate-wise sum of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(std::min(a.size(), b.size()));
		for (size_t i = 0; i < std::min(a.size(), b.size()); i++)
			c.push_back(a[i] + b[i]);
		return c;
	}

	///Coordinate-wise difference of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(std::max(a.size(), b.size()));
		for (size_t i = 0; i < std::min(a.size(), b.size()); i++)
			c.push_back(a[i] - b[i]);
		return c;
	}

	///Coordinate-wise opposite of vector
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a)
	{
		std::vector<T> b;
		b.reserve(a.size());
		for (const auto& i : a)
			b.push_back(-i);
		return b;
	}

	///Product of element and vector
	template <typename T>
	std::vector<T> operator*(T a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(b.size());
		for (const auto& i : b)
			c.push_back(a * i);
		return c;
	}

	///Printing a vector
	template<typename T>
	std::ostream& operator<<(std::ostream& out, const std::vector<T>& a) {
		for (const auto& i : a)
			out << i << ",";
		return out;
	}

	///Least common multiple of vector of elements
	template<typename T>
	int lcm(const T& v) {
		if (v.size() == 0)
			return 0;
		auto n = 1;
		for (const auto& i : v)
			n = std::lcm(n, i);
		return n;
	}


	///Hash vector given minimum and maximum values of its entries.
	template<typename T>
	long hashvector(const T& deg, const T& min, const T& max) {
		long hash = deg[0] - min[0];
		long prod = max[0] - min[0] + 1;
		for (long i = 1; i < deg.size(); i++) {
			hash += (deg[i] - min[i]) * prod;
			prod *= max[i] - min[i] + 1;
		}
		return hash;
	}

	///Turn sparse matrix into a vector of Eigen triplets
	template<typename T, int StorageOrder, typename storage>
	triplets<T, storage> make_triplets(const Eigen::SparseMatrix<T, StorageOrder, storage>& a) {
		triplets<T, storage> b;
		b.reserve(a.nonZeros());
		for (decltype(a.outerSize()) k = 0; k < a.outerSize(); k++) {
			for (typename Eigen::SparseMatrix<T, StorageOrder, storage>::InnerIterator it(a, k); it; ++it)
				b.push_back(Eigen::Triplet<T, storage>(it.row(), it.col(), it.value()));
		}
		return b;
	}


	///Turn sparse matrix into a vector of Eigen triplets
	template<typename T, int StorageOrder, typename storage>
	triplets<T, storage> keep_row_triplets(const Eigen::SparseMatrix<T, StorageOrder, storage>& a, const std::vector<storage>& keep) {
		triplets<T, storage> b;
		b.reserve(a.nonZeros());
		for (decltype(a.outerSize()) k = 0; k < a.outerSize(); k++) {
			for (typename Eigen::SparseMatrix<T, StorageOrder, storage>::InnerIterator it(a, k); it; ++it)
				if (keep[it.row()] != -1)
					b.push_back(Eigen::Triplet<T, storage>(keep[it.row()], it.col(), it.value()));
		}
		return b;
	}
}
