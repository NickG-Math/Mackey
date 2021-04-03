#pragma once
#include "../General.hpp"

///@file
///@brief Contains general functions and vector overloads independent of everything else.

namespace mackey {

	//Constructor using min, max and policy
	template<typename T>
	InterpolatingVectorGenerator<T>::InterpolatingVectorGenerator(const T& min, const T& max) : min(min), max(max) {}

	//Returns the total amount of elements that will be generated
	template<typename T>
	size_t InterpolatingVectorGenerator<T>::size() const {
		if (min.empty())
			return 0;
		size_t total = max[0] - min[0] + 1;
		for (decltype(max.size()) i = 1; i < max.size(); i++)
			total *= max[i] - min[i] + 1;
		return total;
	}

	//Returns vector of all generated elements
	template<typename T>
	std::vector<T> InterpolatingVectorGenerator<T>::operator()() const {
		std::vector<T> v;
		v.reserve(size());
		for (const auto& i : *this)
			v.push_back(i);
		return v;
	}

	//Constant iterator that is used in a ranged for loop to generate the interpolating vectors. Non constant version is illegal
	template<typename T>
	const T& InterpolatingVectorGenerator<T>::const_iterator::operator*() const {
		if (completed == 1){
			std::cerr << "Cannot dereference end() iterator!";
			abort();
		}
		return generated;
	}

	//Inequality of iterators
	template<typename T>
	bool InterpolatingVectorGenerator<T>::const_iterator::operator!=(const const_iterator& other) const {
		if (other.completed) //quick check if other=end
			return !completed;
		else
			return generated != other.generated;
	}

	//Generates next element
	template<typename T>
	auto& InterpolatingVectorGenerator<T>::const_iterator::operator++() {
		decltype(max.size()) i = 0;
		while (i < max.size() && generated[i] == max[i]) {
			generated[i] = min[i];
			i++;
		}
		if (i == max.size())
			completed = 1;
		else
			generated[i]++;
		return *this;
	}

	//Constructor 
	template<typename T>
	InterpolatingVectorGenerator<T>::const_iterator::const_iterator(const T& min, const T& max, bool completed) : min(min), max(max), completed(completed) {
		if (!completed)
			generated = min;
	}

	//Initial generator
	template<typename T>
	auto InterpolatingVectorGenerator<T>::begin() const {
		return const_iterator(min, max, 0);
	}

	//Terminal generator.
	template<typename T>
	auto InterpolatingVectorGenerator<T>::end() const {
		return const_iterator(min, max, 1);
	}

	///Given a vector of vectors v form all combinations of the form {v[0][?],v[1][?],...} for varying ?
	template<typename T>
	std::vector<std::vector<T>> combinations(const std::vector<std::vector<T>>& v) {
		std::vector<std::vector<T>> w;
		std::vector<T> each;
		each.reserve(v.size());
		std::vector<size_t> counter(v.size(), 0);
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
			for (size_t i = 0; i < v.size(); i++)
				each.push_back(v[i][counter[i]]);
			w.push_back(each);
			change = v.size() - 1;
			counter[change]++;
		}
		return w;
	}

	//Given vector or array returns vector starting from index start.
	template<typename T>
	std::vector<T> tail(const T* const& ptr, int size, int start) {
		return std::vector<T>(ptr + start, ptr + size);
	}

	//Makes the multiple of the basis vector 0,...,multiple,...,0
	template<typename rank_t>
	rank_t basisElement(int length, int position, int multiple) {
		rank_t a = rank_t::Zero(length);
		a[position] = multiple;
		return a;
	}

	//Makes a new matrix out of A keeping only its rows indicated by Z
	template<typename T, typename S>
	T KeepRow(const T& A, const S& Z) {
		T B(Z.size(), A.cols());
		size_t j = 0;
		for (const auto& i : Z) {
			B.row(j) = A.row(i);
			j++;
		}
		return B;
	}

	//Makes a new matrix out of A keeping only its columns indicated by Z
	template<typename T, typename S>
	T KeepCol(const T& A, const S& Z) {
		T B(A.rows(), Z.size());
		size_t j = 0;
		for (const auto& i : Z) {
			B.col(j) = A.col(i);
			j++;
		}
		return B;
	}

	//Takes the sum of the entries of an Eigen vector at higher precision than its inputs.
	template<typename Derived>
	int64_t summation(const Eigen::MatrixBase<Derived>& A) {
		return A.template cast<int64_t>().sum();
	}


	template<typename rank_t>
	int64_t summation(const rank_t& u, int64_t limit) {
		int64_t sum = 0;
		for (int64_t i = 0; i < limit; i++)
			sum += u[i];
		return sum;
	}


	//Transforms V[0],V[1],...,V[n] to V[1],V[2],...,V[n],V[0]
	template<typename T>
	void rotate(T& V) {
		auto a = V[0];
		for (int i = 0; i < V.size() - 1; i++)
			V[i] = V[i + 1];
		V[V.size() - 1] = a;
	}


	//Creates a matrix of size m x n that alternates the given vector v
	template<typename Matrix_t>
	Matrix_t altmatrix(int m, int n, const std::vector<scalar_t<Matrix_t>>& v) {
		if constexpr (SFINAE::is_Dense<Matrix_t>::value) {
			Matrix_t matrix(m, n);
			col_vector_t<Matrix_t> alt = Eigen::Map<const col_vector_t<Matrix_t>>(v.data(), v.size());
			col_vector_t<Matrix_t> column(m);
			for (int i = 0; i < m; i += alt.size())
				column.segment(i, alt.size()) = alt;
			for (int j = n - 1; j >= 0; j--) {
				rotate(column);
				matrix.col(j) = column;
			}
			return matrix;
		}
		else {
			auto dense = altmatrix<dense_t<Matrix_t>>(m, n, v);
			return dense.sparseView();
		}
	}



	//Tests if vector is zero
	template <typename T>
	bool isZero(const std::vector<T>& a)
	{
		for (auto i : a) {
			if (i != 0)
				return 0;
		}
		return 1;
	}

	//Coordinate-wise sum of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(std::min(a.size(), b.size()));
		for (size_t i = 0; i < std::min(a.size(), b.size()); i++)
			c.push_back(a[i] + b[i]);
		return c;
	}

	//Coordinate-wise difference of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(std::max(a.size(), b.size()));
		for (size_t i = 0; i < std::min(a.size(), b.size()); i++)
			c.push_back(a[i] - b[i]);
		return c;
	}

	//Coordinate-wise opposite of vector
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a)
	{
		std::vector<T> b;
		b.reserve(a.size());
		for (const auto& i : a)
			b.push_back(-i);
		return b;
	}

	//Product of element and vector
	template <typename T>
	std::vector<T> operator*(T a, const std::vector<T>& b)
	{
		std::vector<T> c;
		c.reserve(b.size());
		for (const auto& i : b)
			c.push_back(a * i);
		return c;
	}


	//Printing a vector
	template<typename T>
	std::ostream& operator<<(std::ostream& out, const std::vector<T>& a) {
		if (a.empty())
			return out;
		for (size_t i = 0; i < a.size() - 1; i++)
			out << a[i] << ",";
		out << a.back();
		return out;
	}
	
	//Least common multiple of vector of elements
	template<typename T>
	int lcm(const T& v) {
		if (v.empty())
			return 0;
		auto n = 1;
		for (auto i : v)
			n = std::lcm(n, i);
		return n;
	}


	//Perfect hash vectors given minimum and maximum values of entries.
	template<typename T, typename S>
	int64_t perfect_hash(const T& deg, const S& min, const S& max) {
		int64_t hash = deg[0] - min[0];
		int64_t prod = 1;
		for (decltype(deg.size()) i = 1; i < deg.size(); i++) {
			prod *= max[i - 1] - min[i - 1] + 1;
			hash += (deg[i] - min[i]) * prod;
		}
		return hash;
	}

	template<typename S>
	size_t Hash<>::operator()(const S& v) const {
		size_t hash = 0;
		for (const auto i : v)
			hash ^= i + 0x9e3779b97f4a7c15 + (hash << 6) + (hash >> 2); //0x9e3779b97f4a7c15 is the golden ratio
		return hash;
	}

	template<typename T>
	Hash<T>::Hash(const T& max) : max(max) {}

	template<typename T>
	size_t Hash<T>::max_hash() const {
		return operator()(max);
	}

	template<typename T>
	template<typename S>
	size_t Hash<T>::operator()(const S& deg) const {
		size_t hash = deg[0];
		size_t prod = 1;
		for (decltype(deg.size()) i = 1; i < deg.size(); i++) {
			prod *= max[i - 1] + 1;
			hash += deg[i] * prod;
		}
		return hash;
	}

	template<typename T>
	Hash<T, T>::Hash(const T& min, const T& max) : min(min), max(max) {}

	template<typename T>
	size_t Hash<T, T>::max_hash() const {
		return operator()(max);
	}
	template<typename T>
	template<typename S>
	size_t Hash<T, T>::operator()(const S& deg) const {
		size_t hash = deg[0] - min[0];
		size_t prod = 1;
		for (decltype(deg.size()) i = 1; i < deg.size(); i++) {
			prod *= max[i - 1] - min[i - 1] + 1;
			hash += (deg[i] - min[i]) * prod;
		}
		return hash;
	}



	//Turn sparse matrix into a vector of Eigen triplets
	template<typename T>
	std::vector<triplet_t<T>> make_triplets(const T& a) {
		std::vector<triplet_t<T>> b;
		b.reserve(a.nonZeros());
		for (decltype(a.outerSize()) k = 0; k < a.outerSize(); k++)
			for (typename T::InnerIterator it(a, k); it; ++it)
				b.emplace_back(it.row(), it.col(), it.value());
		return b;
	}


	//Turn sparse matrix into a vector of Eigen triplets
	template<typename T, typename S>
	std::vector<triplet_t<T>> keep_row_triplets(const T& a, const S& keep) {
		std::vector<triplet_t<T>> b;
		b.reserve(a.nonZeros());
		for (decltype(a.outerSize()) k = 0; k < a.outerSize(); k++) {
			for (typename T::InnerIterator it(a, k); it; ++it)
				if (keep[it.row()] != -1)
					b.emplace_back(keep[it.row()], it.col(), it.value());
		}
		return b;
	}

	namespace implementation_details{

		template<typename T>
		typename T::StorageIndex find_inner_index(typename T::StorageIndex row, typename T::StorageIndex col) {
			if constexpr (SFINAE::is_sparse_row_major<T>::value)
				return col;
			else
				return row;
		}

		template<typename T>
		typename T::StorageIndex find_outer_index(typename T::StorageIndex row, typename T::StorageIndex col) {
			if constexpr (SFINAE::is_sparse_row_major<T>::value)
				return row;
			else
				return col;
		}

		template <typename T>
		class it_nnz_conditional<T, 1>
		{
			typedef typename T::StorageIndex storage;
			typename T::InnerIterator _it;

		protected:
			auto& inner_ind();
			auto& it() { return _it; }
			static constexpr bool is_sparse = 1;
			it_nnz_conditional(const T& matrix, storage row_start, storage col_start)
				: _it(matrix, find_outer_index<T>(row_start, col_start)) {}
		};

		template <typename T>
		class it_nnz_conditional<T, 0>
		{
			typedef typename T::StorageIndex storage;
			storage _inner_ind;

		protected:
			auto& inner_ind() {	return _inner_ind;}
			auto& it();
			static constexpr bool is_sparse = 0;
			it_nnz_conditional(const T&, storage row_start, storage)
				: _inner_ind(row_start) {}
		};
	}

	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	auto IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::value() {
		if constexpr (is_sparse)
			return it().value();
		else
			return matrix(inner_ind(), outer_ind);
	}

	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	auto IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::row() {
		if constexpr (is_sparse)
			return it().row();
		else
			return inner_ind();
	}

	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	auto IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::col() {
		if constexpr (is_sparse)
			return it().col();
		else
			return outer_ind;
	}

	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>& IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::operator++() {
		if constexpr (is_sparse)
			++it();
		else {
			if constexpr (!only_iterate_over_row && !only_iterate_over_col) {
				inner_ind()++;
				if (inner_ind() == matrix.rows()) {
					inner_ind() = inner_start;
					outer_ind++;
				}
			}
			else {
				if constexpr (only_iterate_over_row)
					outer_ind++;
				else
					inner_ind()++;
			}
		}
		increase_until_it_works();
		return *this;
	}

	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::operator bool() {
		if constexpr (is_sparse)
			return it();
		else
			return (inner_ind() < matrix.rows() && outer_ind < matrix.cols());
	}


	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::IteratorNNZ(const T& matrix, storage row_start, storage col_start)
		: implementation_details::it_nnz_conditional<T>(matrix, row_start, col_start),
		matrix(matrix),
		inner_start(implementation_details::find_inner_index<T>(row_start, col_start)),
		outer_ind(implementation_details::find_outer_index<T>(row_start, col_start))
	{
		increase_until_it_works();
	}


	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::IteratorNNZ(const T& matrix, storage start)
		: IteratorNNZ(matrix, start, start) {}


	template<typename T, bool only_iterate_over_row, bool only_iterate_over_col>
	void IteratorNNZ<T, only_iterate_over_row, only_iterate_over_col>::increase_until_it_works() {
		if constexpr (is_sparse) {
			if constexpr (only_iterate_over_row || only_iterate_over_col)
				for (;it() && it().index() < inner_start; ++it()) {}
			else {
				while ((it() && it().index() < inner_start) || (!it() && outer_ind < matrix.outerSize())) {
					for (;it() && it().index() < inner_start; ++it()) {}
					for (;!it() && outer_ind < matrix.outerSize();) {
						outer_ind++;
						if (outer_ind<matrix.outerSize())
							it() = typename T::InnerIterator(matrix, outer_ind);
					}
				}
			}
		}
		else if constexpr (!is_sparse) {
			if constexpr (only_iterate_over_row) {
				for (; outer_ind < matrix.cols(); outer_ind++)
					if (matrix(inner_ind(), outer_ind) != 0)
						break;
			}
			else if constexpr (only_iterate_over_col) {
				for (;inner_ind() < matrix.rows(); inner_ind()++)
					if (matrix(inner_ind(), outer_ind) != 0)
						break;
			}
			else {
				while (outer_ind != matrix.cols() && matrix(inner_ind(), outer_ind) == 0) {
					inner_ind()++;
					if (inner_ind() == matrix.rows()) {
						inner_ind() = inner_start;
						outer_ind++;
					}
				}
			}
		}
	}
}
