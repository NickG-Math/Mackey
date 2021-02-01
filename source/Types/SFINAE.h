#pragma once
#include <iostream>
#include <type_traits>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Coefficients/Z_n.h>

///@file
///@brief Contains all the SFINAE tricks


namespace Mackey {

	//forward declaration for SFINAE
	template<typename T>
	class Graph;

	
	//////////////////
	///SFINAE is used to detect types and existence of methods
	
	///The test functions shouldn't be public, but Clang complains otherwise for some reason.
	//////////////////
	struct SFINAE {

		//Detect Eigen Dense matrices 
		static constexpr std::false_type test_Dense(...);

		template<typename T>
		static constexpr std::true_type test_Dense(Eigen::MatrixBase<T>*);

		//Detect Eigen Sparse matrices
		static constexpr	std::false_type test_Sparse(...);

		template <typename T>
		static constexpr std::true_type test_Sparse(Eigen::SparseMatrixBase<T>*);

		//Detect Eigen Sparse Row Major matrices
		static constexpr	std::false_type test_SRM(...);

		template <typename T>
		static constexpr std::true_type test_SRM(Eigen::SparseMatrix<T,1>*);

		///Tests if matrix is dense
		template<typename T>
		using is_Dense = decltype(test_Dense(std::declval<T*>()));

		///Tests if matrix is sparse
		template<typename T>
		using is_Sparse = decltype(test_Sparse(std::declval<T*>()));

		///Tests if matrix is sparse row major
		template<typename T>
		using is_sparse_row_major = decltype(test_SRM(std::declval<T*>()));


		//Detect if T=Z<N>
		static constexpr std::false_type test_finite_cyclic(...);

		template <int N>
		static constexpr std::true_type test_finite_cyclic(Z<N>);

		///Tests if T=Z<N>
		template<typename T>
		using is_finite_cyclic = decltype(test_finite_cyclic(std::declval<T>()));


		//Detect implementation of computePath()
		template<typename T>
		static constexpr std::false_type test_computePath(...);

		template<typename T>
		static constexpr decltype(std::declval<T>().computePath(), std::true_type()) test_computePath(int);

		///Tests if T has member computePath()
		template<typename T>
		using has_computePath = decltype(test_computePath<T>(0));


		//Detect implementation of initialize()

		template<typename T>
		static constexpr std::false_type test_initialize(...);

		template<typename T>
		static constexpr decltype(std::declval<T>().initialize(), std::true_type()) test_initialize(int);

		///Tests if T has member initialize()
		template<typename T>
		using has_initialize = decltype(test_initialize<T>(0));


		///Detects graphs
		static constexpr std::false_type is_graph(...);

		template <typename T>
		static constexpr std::true_type is_graph(Graph<T>);

		///Tests if T is a graph type
		template <typename T>
		using is_graph_t = decltype(is_graph(std::declval<T>()));

	};


}
