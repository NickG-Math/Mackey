#pragma once
#include <vector>
#include "Types/Aliases.hpp"
#include <numeric> //lcm
#include <iostream> //cerr

///@file
///@brief Contains general functions and vector overloads independent of everything else.

///Everything in this library is under this namespace
namespace mackey
{

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief		Generates all vectors interpolating between given min and max vectors
	///	@tparam	T	The type of vectors eg \c std::vector<int>
	///	@details	Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a vector_interpolate_generator object.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename T>
	class InterpolatingVectorGenerator
	{
	public:
		///	@brief	Constructor using min, max and policy
		InterpolatingVectorGenerator(const T& min, const T& max);

		///	@brief	Returns the total amount of elements that will be generated
		size_t size() const;

		///	@brief	Returns vector of all generated elements
		std::vector<T> operator()() const;

		///	@brief	Constant iterator through generated vectors
		///	@details Used in a ranged for loop to generate the interpolating vectors. Non constant version is illegal
		class const_iterator
		{
		public:
			///Returns the generated vector
			const T& operator*() const;
			///Inequality of iterators
			bool operator!=(const const_iterator& other) const;
			///Generates next element
			auto& operator++();

		private:
			const_iterator(const T& min, const T& max, bool completed);
			const T& min, max;
			bool completed;
			T generated;
			friend class InterpolatingVectorGenerator;
		};
		///	@brief	Initial generator
		auto begin() const;
		///	@brief	Terminal generator.
		auto end() const;

	private:
		const T& min, max;
	};

	/// @brief		Given a vector of vectors \c v form all combinations of the form {v[0][?],v[1][?],...} for varying ?
	/// @tparam T	Value type of vector of vectors \c v
	/// @param v	The vector of vector v
	/// @return		A vector of all vectors that are of the form {v[0][?],v[1][?],...} for varying ?
	template <typename T>
	std::vector<std::vector<T>> combinations(const std::vector<std::vector<T>>& v);

	///	@brief	Given vector or array returns vector starting from index start.
	///	@details	Example: Given v apply as tail(v.data(),v.size(), start)
	template <typename T>
	std::vector<T> tail(const T* const& ptr, int size, int start);

	///Makes the multiple of the basis vector 0,...,multiple,...,0
	template <typename rank_t>
	rank_t basisElement(int length, int position, int multiple = 1);

	///Makes a new matrix out of A keeping only its rows indicated by Z
	template <typename T, typename S>
	T KeepRow(const T& A, const S& Z);

	///Makes a new matrix out of A keeping only its columns indicated by Z
	template <typename T, typename S>
	T KeepCol(const T& A, const S& Z);

	///Takes the sum of the entries of an Eigen vector at higher precision than its inputs.
	template <typename Derived>
	int64_t summation(const Eigen::MatrixBase<Derived>& A);

	template <typename rank_t>
	int64_t summation(const rank_t& u, int64_t limit);

	///Transforms V[0],V[1],...,V[n] to V[1],V[2],...,V[n],V[0]
	template <typename T>
	void rotate(T& V);

	///	@brief	Creates a matrix of size m x n that alternates the given vector v
	///	@details Eg for m=2, n=2, v={a,b} we get {a,b // b,a}
	template <typename Matrix_t>
	Matrix_t altmatrix(int m, int n, const std::vector<scalar_t<Matrix_t>>& v);

	///Tests if vector is zero
	template <typename T>
	bool isZero(const std::vector<T>& a);

	///Coordinate-wise sum of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b);

	///Coordinate-wise difference of vectors (up to the minimum of their lengths)
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b);

	///Coordinate-wise opposite of vector
	template <typename T>
	std::vector<T> operator-(const std::vector<T>& a);

	///Product of element and vector
	template <typename T>
	std::vector<T> operator*(T a, const std::vector<T>& b);

	///Printing a vector
	template <typename T>
	std::ostream& operator<<(std::ostream& out, const std::vector<T>& a);

	///Least common multiple of vector of elements
	template <typename T>
	int lcm(const T& v);

	///General template class for hashing vectors 
	template <typename...>
	class Hash;

	///Hashes vectors with collisions
	template <>
	class Hash<>
	{
	public:
		template <typename S>
		size_t operator()(const S&) const;
	};

	///Perfect hashes vectors between 0 and given max
	template <typename T>
	class Hash<T>
	{
		const T max;

	public:
		Hash(const T&);
		size_t max_hash() const;
		template <typename S>
		size_t operator()(const S&) const;
	};

	///Perfect hashes vectors between given min and max
	template <typename T>
	class Hash<T, T>
	{
		const T min, max;

	public:
		size_t max_hash() const;
		Hash(const T&, const T&);
		template <typename S>
		size_t operator()(const S&) const;
	};

	///Perfect hash vectors given minimum and maximum values of entries.
	template <typename T, typename S>
	int64_t perfect_hash(const T& deg, const S& min, const S& max);

	///Turn sparse matrix into a vector of Eigen triplets
	template <typename T>
	std::vector<triplet_t<T>> make_triplets(const T& a);

	///Turn sparse matrix into a vector of Eigen triplets
	template <typename T, typename S>
	std::vector<triplet_t<T>> keep_row_triplets(const T& a, const S& keep);

	namespace implementation_details
	{
		///	@brief Helper class for iterating over nonzero entries in a sparse or dense Eigen matrix (offering consistent API) 
		template <typename T, bool = mackey::SFINAE::is_Sparse<T>::value>
		class it_nnz_conditional;
	}

	//////////////////////////////////////////////////////////////////////////////
	///	@brief				Iterator over nonzero entries of Eigen matrices
	///	@tparam	T			The type of the Eigen matrix
	///	@tparam	row_only	Set to 1 to iterate only through the entries in the given row (so the block is a row)
	///	@tparam	col_only	Set to 1 to iterate only through the entries in the given column (so the block is a column)
	///	@note				If both boolean templates are 0 then the block is any lower right corner
	/////////////////////////////////////////////////////////////////////////////
	template <typename T, bool row_only, bool col_only>
	class IteratorNNZ : implementation_details::it_nnz_conditional<T>
	{
		typedef typename T::StorageIndex storage;
		using implementation_details::it_nnz_conditional<T>::is_sparse;
		using implementation_details::it_nnz_conditional<T>::inner_ind;
		using implementation_details::it_nnz_conditional<T>::it;
	public:
		auto value();							///<The entry of the matrix
		auto row();								///<The row index
		auto col();								///<The column index
		IteratorNNZ& operator++(); ///<Increase iterator
		operator bool();						///<Is 1 if iterator is valid, 0 if nonvalid

		/// @brief				Constructor given matrix and block (lower right corner)
		/// @param matrix		The matrix to iterate through
		/// @param row_start	The starting row of the block: we only visit nonzero entries \c matrix(i,j) with i>=row_start
		/// @param col_start	The starting column of the block: we only visit nonzero entries \c matrix(i,j) with j>=col_start
		IteratorNNZ(const T& matrix, storage row_start, storage col_start);

		/// @brief			Constructor given matrix and symmetric block
		/// @param matrix	The matrix to iterate through
		/// @param start	The starting row+column of the block: we only visit nonzero entries \c matrix(i,j) with i,j>=start
		IteratorNNZ(const T& matrix, storage start);

	private:
		const T& matrix;
		const storage inner_start;
		storage outer_ind;
		void increase_until_it_works();
	};
}
#include "impl/General.ipp"