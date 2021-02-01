#pragma once

///Enables serialization of certain private members
#define CEREALIZE

#include <fstream>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp> //std::pair
#include <cereal/types/string.hpp>
#include "Additive.h"
#include "Factorization.h"

///@file
///@brief Contains the methods for serializing the results of our computations (before their interpretation) via cereal

///Consult the Cereal documentation for this
namespace cereal

{	///Save Eigen dense matrix
	template<typename Archive, typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
	void save(Archive& archive, const Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>& A)
	{
		long rows = A.rows();
		long cols = A.cols(); //useful to have both if A has one 0 dimension
		std::vector<T> vec(A.data(), A.data() + A.size());
		archive(CEREAL_NVP(rows), CEREAL_NVP(cols), CEREAL_NVP(vec));
	}


	///Load Eigen dense matrix
	template<typename Archive, typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
	void load(Archive& archive, Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>& A)
	{
		long rows, cols;
		std::vector<T> vec;
		archive(CEREAL_NVP(rows), CEREAL_NVP(cols), CEREAL_NVP(vec));
		A = Eigen::Map<Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>>(vec.data(), rows, cols);
	}


	/////Serialize triplets
	//template<typename Archive, typename T>
	//void save(Archive& archive, const Eigen::Triplet<T>& A)
	//{
	//	int row = A.row();
	//	int col = A.col();
	//	T value = A.value();
	//	archive(CEREAL_NVP(row), CEREAL_NVP(col), CEREAL_NVP(value));
	//}


	/////Serialize triplets
	//template<typename Archive, typename T>
	//void load(Archive& archive, Eigen::Triplet<T>& A)
	//{
	//	int row, col;
	//	T value;
	//	archive(CEREAL_NVP(row), CEREAL_NVP(col), CEREAL_NVP(value));
	//	A = Eigen::Triplet<T>(row, col, value);
	//}



	/////Save Eigen sparse matrix
	//template<typename Archive, typename T, int StorageOrder>
	//void save(Archive& archive, const Eigen::SparseMatrix<T, StorageOrder>& A)
	//{
	//	triplets<T> B = make_triplets(A);
	//	archive(CEREAL_NVP(B));
	//}


	/////Load Eigen sparse matrix
	//template<typename Archive, typename T, int StorageOrder>
	//void load(Archive& archive, Eigen::SparseMatrix<T, StorageOrder>& A)
	//{
	//	triplets<T> B;
	//	archive(CEREAL_NVP(B));
	//	A.setFromTriplets(B.begin(), B.end());
	//}

}
namespace Mackey {

	///Chains cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, Chains<rank_t, diff_t>& C) {
		archive(CEREAL_NVP(C.maxindex), CEREAL_NVP(C.rank), CEREAL_NVP(C.diff));
	}

	///Chains cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, Homology<rank_t, diff_t>& H) {
		archive(CEREAL_NVP(H.Groups), CEREAL_NVP(H.Generators), CEREAL_NVP(H.isZero), CEREAL_NVP(H.nonZeroVectors), CEREAL_NVP(H.dontModOut), \
			CEREAL_NVP(H.In_Q), CEREAL_NVP(H.Out_Qi), CEREAL_NVP(H.In_P_full), CEREAL_NVP(H.In_P_reduced));
	}

	///IDGenerator cerealize
	template<typename Archive, typename rank_t>
	void serialize(Archive& archive, IDGenerators<rank_t>& ID) {
		archive(CEREAL_NVP(ID.group), CEREAL_NVP(ID.group_lower), CEREAL_NVP(ID.Tr), CEREAL_NVP(ID.Res));
	}

	///MackeyFunctor cerealize
	template<typename Archive, typename rank_t>
	void serialize(Archive& archive, MackeyFunctor<rank_t>& Mack) {
		archive(CEREAL_NVP(Mack.Groups), CEREAL_NVP(Mack.Tr), CEREAL_NVP(Mack.Res), CEREAL_NVP(Mack.Weyl), CEREAL_NVP(Mack.name));
	}


	///AdditiveStructure cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, AdditiveStructure<rank_t, diff_t>& A) {
		archive(CEREAL_NVP(A.minsphere), CEREAL_NVP(A.maxsphere), CEREAL_NVP(A.allMackeys));
	}

	///Green cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, Green<rank_t, diff_t>& G) {
		archive(CEREAL_NVP(G.Groups), CEREAL_NVP(G.basis), CEREAL_NVP(G.boxID), CEREAL_NVP(G.first_number_selections), CEREAL_NVP(G.second_number_selections), CEREAL_NVP(G.isZero));
	}

	///MultiplicationTable cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void save(Archive& archive, const MultiplicationTable<rank_t, diff_t>& M) {
		archive(CEREAL_NVP(M.level), CEREAL_NVP(M.NonZeroHomology), CEREAL_NVP(M.degree), CEREAL_NVP(M.antidegree), CEREAL_NVP(M.index_product), \
			CEREAL_NVP(M.minsphere), CEREAL_NVP(M.maxsphere), CEREAL_NVP(M.Greens), CEREAL_NVP(M.basicIrreducibles), CEREAL_NVP(M.number_of_irreducibles), CEREAL_NVP(M.tripleGreens));
	}

	///MultiplicationTable cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void load(Archive& archive, MultiplicationTable<rank_t, diff_t>& M) {
		archive(CEREAL_NVP(M.level), CEREAL_NVP(M.NonZeroHomology), CEREAL_NVP(M.degree), CEREAL_NVP(M.antidegree), CEREAL_NVP(M.index_product), \
			CEREAL_NVP(M.minsphere), CEREAL_NVP(M.maxsphere), CEREAL_NVP(M.Greens), CEREAL_NVP(M.basicIrreducibles), CEREAL_NVP(M.number_of_irreducibles), CEREAL_NVP(M.tripleGreens));
		M.getIrreducibleChains();
	}


	///Saves object to binary file of given name. Serialization is provided by cereal
	template<typename T>
	void saver(const T& A, const std::string& filename, const std::string& type) {
		std::cout << "Saving to " << filename << "\n";
		if (type == "binary") {
			std::ofstream ofs(filename, std::ofstream::binary);
			cereal::BinaryOutputArchive oarchive(ofs);
			oarchive(A);
		}
		else if (type == "xml") {
			std::ofstream ofs(filename);
			cereal::XMLOutputArchive oarchive(ofs);
			oarchive(A);
		}
		else if (type == "json") {
			std::ofstream ofs(filename);
			cereal::JSONOutputArchive oarchive(ofs);
			oarchive(A);
		}
		std::cout << "Saved \n";
	}

	///Loads object from binary file of given name. Serialization is provided by cereal
	template<typename T>
	void loader(T& A, const std::string& filename, const std::string& type) {
		std::cout << "Loading from " << filename << "\n";
		if (type == "binary") {
			std::ifstream ifs(filename, std::ifstream::binary);
			cereal::BinaryInputArchive iarchive(ifs);
			iarchive(A);
		}
		else {
			std::ifstream ifs(filename);
			if (type == "xml") {
				cereal::XMLInputArchive iarchive(ifs);
				iarchive(A);
			}
			else if (type == "json") {
				cereal::JSONInputArchive iarchive(ifs);
				iarchive(A);
			}
		}
		std::cout << "Loaded \n";
	}
}
