#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "SFINAE.hpp"

///@file
///@brief Contains general template aliases. 


namespace mackey {

	///Scalar of matrix
	template<typename T>
	using scalar_t = typename T::Scalar;

	///Storage type of matrix
	template<typename T>
	using storage_t = typename T::StorageIndex;

	///Dense row matrix
	template<typename T>
	using row_vector_t = Eigen::Matrix<scalar_t<T>,1,-1>;
	
	///Dense column matrix
	template<typename T>
	using col_vector_t = Eigen::Matrix<scalar_t<T>, -1, 1>;

	///Dense column major
	template<typename T>
	using dense_t= Eigen::Matrix<scalar_t<T>,-1,-1>;

	///Sparse column major
	template<typename T>
	using sparse_t = Eigen::SparseMatrix<scalar_t<T>, 0, storage_t<T>>;

	template<typename T>
	using row_major_t = std::conditional<SFINAE::is_Dense<T>::value, Eigen::Matrix<scalar_t<T>, -1, -1,1>, Eigen::SparseMatrix<scalar_t<T>, 1, storage_t<T>>>;

	template<typename T>
	using col_major_t = std::conditional<SFINAE::is_Dense<T>::value, Eigen::Matrix<scalar_t<T>, -1, -1>, Eigen::SparseMatrix<scalar_t<T>, 0, storage_t<T>>>;

	///Eigen triplets type
	template<typename T>
	using triplet_t = Eigen::Triplet<scalar_t<T>, storage_t<T>>;

	template<typename, typename>
	struct Arrow;

	template<typename, typename>
	class Chains;

	template<typename, typename>
	class Junction;

	///Arrow given group_t
	template<typename T>
	using arrow_t = Arrow<typename T::rank_t, typename T::diff_t>;

	///Chains given group_t
	template<typename T>
	using chains_t = Chains<typename T::rank_t, typename T::diff_t>;

	///Junction given group_t
	template<typename T>
	using junction_t = Junction<typename T::rank_t, typename T::diff_t>;

	template<typename, typename>
	class Homology;

	///Type of generators in homology
	template<typename rank_t, typename diff_t>
	using gen_t = typename Homology<rank_t, diff_t>::gen_t;

}
