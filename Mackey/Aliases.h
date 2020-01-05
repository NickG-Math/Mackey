#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

///@file
///@brief Contains general template aliases. 


namespace Mackey {

	///Scalar of matrix
	template<typename T>
	using Scalar_t = typename T::Scalar;

	///Row matrix
	template<typename T>
	using row_t = typename Eigen::Matrix<Scalar_t<T>,1,-1>;
	
	///Column matrix
	template<typename T>
	using col_t = typename Eigen::Matrix<Scalar_t<T>, -1, 1>;

	///Column major
	template<typename T>
	using mat_t= typename Eigen::Matrix<Scalar_t<T>,-1,-1>;

	///Row major
	template<typename T>
	using mat_t_r = typename Eigen::Matrix<Scalar_t<T>, -1, -1, Eigen::RowMajor>;

	///Column major sparse
	template<typename T>
	using spm_t = typename Eigen::SparseMatrix<Scalar_t<T>, 0>;

	///Row major sparse
	template<typename T>
	using spm_t_r = typename Eigen::SparseMatrix<Scalar_t<T>, 1>;

	template<typename, typename>
	class Homology;


	///Type of Generator matrices in homology
	template<typename rank_t, typename diff_t>
	using Gen_t = typename Homology<rank_t, diff_t>::Gen_t;

	///Type of generators in homology
	template<typename rank_t, typename diff_t>
	using gen_t = typename Homology<rank_t, diff_t>::gen_t;

}