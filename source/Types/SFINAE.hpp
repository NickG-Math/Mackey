#pragma once
#include <type_traits>
///	@file
///	@brief Contains all the SFINAE tricks


namespace Eigen {
	//forward declaration for SFINAE
	template<typename>
	class MatrixBase;

	template<typename>
	class SparseMatrixBase;

	template<typename, int, typename>
	class SparseMatrix;
}


namespace mackey {

	//forward declaration for SFINAE
	template<int64_t, typename>
	class Z_mod;

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief		Contains typedefs used to detect types and existence of methods
	///	@todo		Replace with concepts
	///	@details	Only use the typedefs eg  constexpr bool b= SFINAE::is_Sparse<T>::value;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	namespace SFINAE {

		//If the test functions are not public then Graph CRTP may fail leading to segfault
	
		///Detects if T is a dense matrix
		constexpr std::false_type test_Dense(...);

		///Detects if T is a dense matrix
		template<typename T>
		constexpr std::true_type test_Dense(Eigen::MatrixBase<T>*);

		///Detects if T is a sparse matrix
		constexpr std::false_type test_Sparse(...);

		///Detects if T is a sparse matrix
		template <typename T>
		constexpr std::true_type test_Sparse(Eigen::SparseMatrixBase<T>*);

		///Detects if T is sparse row major matrix
		constexpr std::false_type test_SRM(...);

		///Detects if T is sparse row major matrix
		template <typename T, typename S>
		constexpr std::true_type test_SRM(Eigen::SparseMatrix<T,1,S>*);

		///Detects if T=Z<N>
		constexpr std::false_type test_finite_cyclic(...);

		///Detects if T=Z<N>
		template <int64_t N, typename T>
		constexpr std::true_type test_finite_cyclic(Z_mod<N,T>);

		///Detects if T has an operator()(int,int)
		template<typename T>
		constexpr std::false_type test_operator_evaluation_int_int(...);

		///Detects if T has an operator()(int,int)
		template<typename T>
		constexpr decltype(std::declval<T>().operator()(std::declval<int>(), std::declval<int>()), std::true_type()) test_operator_evaluation_int_int(int);
	
		///Detects if T is a dense matrix
		template<typename T>
		using is_Dense = decltype(test_Dense(std::declval<T*>()));

		///Detects if T is a sparse matrix
		template<typename T>
		using is_Sparse = decltype(test_Sparse(std::declval<T*>()));

		///Detects if T is sparse row major matrix
		template<typename T>
		using is_sparse_row_major = decltype(test_SRM(std::declval<T*>()));

		///Detects if T=Z<N>
		template<typename T>
		using is_finite_cyclic = decltype(test_finite_cyclic(std::declval<T>()));

		///Detects if T has an operator()(int,int)
		template<typename T>
		using has_operator_evaluation_int_int = decltype(test_operator_evaluation_int_int<T>(0));

	};


}
