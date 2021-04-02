#pragma once

#include <fstream>
#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp> //std::pair
#include <cereal/types/string.hpp>
#include "Coefficients/Z_n.hpp"
#include "Mackey_Functors/Additive.hpp"
#include "Factorization/Factorization.hpp"

///@file
///@brief Contains the methods for serializing the results of our computations (before their interpretation) via cereal

///Consult the Cereal documentation for this
namespace cereal
{	///Save Eigen dense matrix
	template<typename Archive, typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
	void save(Archive& archive, const Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>& A);

	///Load Eigen dense matrix
	template<typename Archive, typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
	void load(Archive& archive, Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>& A);

	///Save Eigen sparse matrix
	template<typename Archive, typename T, int StorageOrder, typename storage_t>
	void save(Archive& archive, const Eigen::SparseMatrix<T, StorageOrder, storage_t>& A);


	///Load Eigen sparse matrix
	template<typename Archive, typename T, int StorageOrder, typename storage_t>
	void load(Archive& archive, Eigen::SparseMatrix<T, StorageOrder, storage_t>& A);

}
namespace mackey {

	///Serialize Z2
	template<typename Archive>
	void serialize(Archive& archive, Z2& A);

	///Chains cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, Chains<rank_t, diff_t>& C);

	///Chains cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, Homology<rank_t, diff_t>& H);

	///IDGenerator cerealize
	template<typename Archive, typename rank_t>
	void serialize(Archive& archive, IDGenerators<rank_t>& ID);

	///MackeyFunctor cerealize
	template<typename Archive, typename rank_t>
	void serialize(Archive& archive, MackeyFunctor<rank_t>& Mack);

	///AdditiveStructure cerealize
	template<typename Archive, typename group_t>
	void serialize(Archive& archive, AdditiveStructure<group_t>& A);

	///Green cerealize
	template<typename Archive, typename group_t>
	void serialize(Archive& archive, Green<group_t>& G);

	///Mult_Table cerealize
	template<typename Archive, typename group_t>
	void serialize(Archive& archive, MultTableData<group_t>& M);

	///Saves object to binary file of given name. Serialization is provided by cereal
	template<typename T>
	void save(const T& A, const std::string& filename);

	///Loads object from binary file of given name. Serialization is provided by cereal
	template<typename T>
	void load(T& A, const std::string& filename);
}
#include "impl/Cerealizer.ipp"
