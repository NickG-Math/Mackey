#pragma once
#include "../Cerealizer.hpp"

///@file
///@brief Contains the methods for serializing the results of our computations (before their interpretation) via cereal

///Consult the Cereal documentation for this
namespace cereal

{	///Save Eigen dense matrix
	template<typename Archive, typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
	void save(Archive& archive, const Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>& A)
	{
		int64_t rows = A.rows();
		int64_t cols = A.cols(); //useful to have both if A has one 0 dimension
		std::vector<T> vec(A.data(), A.data() + A.size());
		archive(CEREAL_NVP(rows), CEREAL_NVP(cols), CEREAL_NVP(vec));
	}


	///Load Eigen dense matrix
	template<typename Archive, typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
	void load(Archive& archive, Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>& A)
	{
		int64_t rows, cols;
		std::vector<T> vec;
		archive(CEREAL_NVP(rows), CEREAL_NVP(cols), CEREAL_NVP(vec));
		A = Eigen::Map<Eigen::Matrix<T, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>>(vec.data(), rows, cols);
	}



	///Serialize triplets
	template<typename Archive, typename T>
	void load(Archive& archive, Eigen::Triplet<T>& A)
	{
		int64_t row, col;
		T value;
		archive(CEREAL_NVP(row), CEREAL_NVP(col), CEREAL_NVP(value));
		A = Eigen::Triplet<T>(row, col, value);
	}

	///Serialize triplets
	template<typename Archive, typename T>
	void save(Archive& archive, const Eigen::Triplet<T>& A)
	{
		int64_t row = A.row();
		int64_t col = A.col();
		T value = A.value();
		archive(CEREAL_NVP(row), CEREAL_NVP(col), CEREAL_NVP(value));
	}



	///Save Eigen sparse matrix
	template<typename Archive, typename T, int StorageOrder, typename storage_t>
	void save(Archive& archive, const Eigen::SparseMatrix<T, StorageOrder, storage_t>& A)
	{
		std::vector<Eigen::Triplet<T, storage_t>> B;
		storage_t rows, cols;
		B.reserve(A.nonZeros());
		for (mackey::IteratorNNZ< Eigen::SparseMatrix<T, StorageOrder, storage_t>, 0, 0> it(A, 0); it; ++it)
			B.emplace_back(it.row(), it.col(), it.value());
		rows = A.rows();
		cols = A.cols();
		archive(CEREAL_NVP(rows), CEREAL_NVP(cols), CEREAL_NVP(B));
	}


	///Load Eigen sparse matrix
	template<typename Archive, typename T, int StorageOrder, typename storage_t>
	void load(Archive& archive, Eigen::SparseMatrix<T, StorageOrder, storage_t>& A)
	{
		std::vector<Eigen::Triplet<T, storage_t>> B;
		storage_t rows, cols;
		archive(CEREAL_NVP(rows), CEREAL_NVP(cols), CEREAL_NVP(B));
		A.resize(rows, cols);
		A.setFromTriplets(B.begin(), B.end());
	}

}
namespace mackey {

	///Serialize Z2
	template<typename Archive>
	void serialize(Archive& archive, Z2& A)
	{
		archive(CEREAL_NVP(A.x));
	}

	///Chains cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, Chains<rank_t, diff_t>& C) {
		archive(CEREAL_NVP(C.rank), CEREAL_NVP(C.diff));
	}

	///Chains cerealize
	template<typename Archive, typename rank_t, typename diff_t>
	void serialize(Archive& archive, Homology<rank_t, diff_t>& H) {
		archive(CEREAL_NVP(H.group), CEREAL_NVP(H.generators), CEREAL_NVP(H.isZero), CEREAL_NVP(H.nonZeroVectors), CEREAL_NVP(H.dontModOut), \
			CEREAL_NVP(H.In_Q), CEREAL_NVP(H.Out_Qi), CEREAL_NVP(H.In_P_full), CEREAL_NVP(H.In_P_reduced));
	}

	///AbelianGroup cerealize
	template<typename Archive, typename T>
	void serialize(Archive& archive, AbelianGroup<T>& Ab) {
		archive(CEREAL_NVP(Ab.group));
	}


	///IDGenerator cerealize
	template<typename Archive, typename rank_t>
	void serialize(Archive& archive, IDGenerators<rank_t>& ID) {
		archive(CEREAL_NVP(ID.group), CEREAL_NVP(ID.group_lower), CEREAL_NVP(ID.tr), CEREAL_NVP(ID.res));
	}

	///MackeyFunctor cerealize
	template<typename Archive, typename rank_t>
	void serialize(Archive& archive, MackeyFunctor<rank_t>& Mack) {
		archive(CEREAL_NVP(Mack.group), CEREAL_NVP(Mack.tr), CEREAL_NVP(Mack.res), CEREAL_NVP(Mack.act), CEREAL_NVP(Mack.name));
	}


	///AdditiveStructure cerealize
	template<typename Archive, typename group_t>
	void serialize(Archive& archive, AdditiveStructure<group_t>& A) {
		archive(CEREAL_NVP(A.minsphere), CEREAL_NVP(A.maxsphere), CEREAL_NVP(A.allMackeys));
	}

	///Green cerealize
	template<typename Archive, typename group_t>
	void serialize(Archive& archive, Green<group_t>& G) {
		archive(CEREAL_NVP(G.group), CEREAL_NVP(G.basis), CEREAL_NVP(G.boxID), CEREAL_NVP(G.first_number_selections), CEREAL_NVP(G.second_number_selections), CEREAL_NVP(G.isZero));
	}

	///Mult_Table cerealize
	template<typename Archive, typename group_t>
	void serialize(Archive& archive, MultTableData<group_t>& M) {
		archive
		(CEREAL_NVP(M.level), CEREAL_NVP(M.NonZeroHomology), CEREAL_NVP(M.degree), CEREAL_NVP(M.antidegree), CEREAL_NVP(M.index_product),
			CEREAL_NVP(M.minsphere), CEREAL_NVP(M.maxsphere), CEREAL_NVP(M.Greens), CEREAL_NVP(M.basicIrreducibles), CEREAL_NVP(M.tripleGreens)
		);
	}


	///Saves object to binary file of given name. Serialization is provided by cereal
	template<typename T>
	void save(const T& A, const std::string& filename) {
		std::cout << "Saving to " << filename << "\n";
		std::ofstream ofs(filename, std::ofstream::binary);
		cereal::BinaryOutputArchive oarchive(ofs);
		oarchive(A);
		std::cout << "Saved \n";
	}

	///Loads object from binary file of given name. Serialization is provided by cereal
	template<typename T>
	void load(T& A, const std::string& filename) {
		std::cout << "Loading from " << filename << "\n";
		std::ifstream ifs(filename, std::ifstream::binary);
		cereal::BinaryInputArchive iarchive(ifs);
		iarchive(A);
		std::cout << "Loaded \n";
	}
}
