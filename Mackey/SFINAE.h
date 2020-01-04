#pragma once
#include <iostream>
#include <type_traits>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Z_n.h"

///@file
///@brief Contains all the SFINAE tricks


namespace {

	//Detect Eigen Dense matrices 
	std::false_type test_Dense(...) { return std::false_type(); }
	template <typename T>
	std::true_type test_Dense(Eigen::MatrixBase<T>) { return std::true_type(); }

	//Detect Eigen Sparse matrices
	std::false_type test_Sparse(...) { return std::false_type(); }
	template <typename T>
	std::true_type test_Sparse(Eigen::SparseMatrixBase<T>) { return std::true_type(); }


	//Detect if T=Z<N>
	std::false_type test_finite_cyclic(...) { return std::false_type(); }
	template <int N>
	std::true_type test_finite_cyclic(Z<N>) { return std::true_type(); }

}

namespace Mackey {

	///Detect Eigen Dense matrices
	template<typename T>
	struct is_Dense : public decltype(test_Dense(std::declval<T>())) {};

	///Detect Eigen Sparse matrices
	template<typename T>
	struct is_Sparse : public decltype(test_Sparse(std::declval<T>())) {};

	///Detects finite coefficients Z<N>
	template<typename T>
	struct is_Finite_Cyclic : public decltype(test_finite_cyclic(std::declval<T>())) {};


}
