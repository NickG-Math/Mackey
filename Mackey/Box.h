#pragma once
#include "General.h" 
#include "Chains.h" 
#include "ChangeBasis.h"

///@file
///@brief Contains the construction of box products of Chains.


namespace Mackey {
	//////////////////////////////////////////////
	/// The box (tensor) product of Chains at a single index

	/// We only store a single differential and the ranks of its domain and range
	//////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	class BoxPoint {
	public:
		rank_t rank_domain;		///<The rank of the domain of the differential.
		rank_t rank_range;	///<The rank of the range of the differential.
		diff_t diff; ///<The differential.

		/// The detailed rank retains the layout of the tensor product. That's used for forming products of generators (padding)
		std::vector<long> detailedrank;

		///Construct from chains C,D the boxpoint (C box D)_i->(C box D)_i-1
		BoxPoint(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int i)
			: C(C), D(D), i(i), lowlimit(std::max(0, i - D.maxindex)), highlimit(std::min(i, C.maxindex)) {
			getranks();
			if (rank_domain.size() != 0 && rank_range.size() != 0)
				getdiff();
		}

	private:
		const Chains<rank_t, diff_t>& C, D;
		const int i, lowlimit, highlimit;
		void getranks();
		void getdiff();
	};


	///The tensor product of ranks of modules A,B
	template<typename rank_t>
	rank_t rankmult(const rank_t& A, const rank_t& B) {
		int sum = 0;
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < B.size(); j++) {
				sum += std::min(A(i), B(j));
			}
		}
		rank_t C(sum);
		int current = 0;
		for (int j = 0; j < B.size(); j++) {
			for (int i = 0; i < A.size(); i++) {
				C.segment(current, std::min(A(i), B(j))).setConstant(std::max(A(i), B(j)));
				current += std::min(A(i), B(j));
			}
		}
		return C;
	}

	///The rank of (C box D)_i
	template<typename rank_t>
	std::pair<rank_t, std::vector<long>> rankBox(const std::vector<rank_t>& C, const std::vector<rank_t>& D, int i) {
		rank_t rank;
		std::vector<long> detailedrank;
		std::vector<rank_t> verydetailedrank;
		detailedrank.reserve(i + 1);
		verydetailedrank.resize(i + 1);
		auto lowlimit = std::max(0, i - (int)D.size()+1);
		auto highlimit = std::min(i, (int)C.size()-1);
		long sum = 0;
		for (int j = lowlimit; j <= highlimit; j++) {
			verydetailedrank[j] = rankmult(C[j], D[i - j]);
			sum += verydetailedrank[j].size();
		}
		rank.resize(sum);
		long track = 0;
		for (const auto& i : verydetailedrank) {
			rank.segment(track, i.size()) = i;
			detailedrank.push_back(summation(i));
			track += i.size();
		}
		return std::make_pair(rank, detailedrank);
	}


	///The rank of (C box D)
	template<typename rank_t>
	std::vector<rank_t> rankBox(const std::vector<rank_t>& C, const std::vector<rank_t>& D) {
		std::vector<rank_t> rank;
		rank.reserve(C.size() + D.size() - 1);
		for (int i = 0; i<C.size()+D.size()-1; i++)
			rank.push_back(rankBox(C,D,i).first);
		return rank;
	}


	template<typename rank_t, typename diff_t>
	void BoxPoint<rank_t, diff_t> ::getranks() {
		auto Ranks = rankBox(C.rank, D.rank, i);
		rank_domain = Ranks.first;
		detailedrank = Ranks.second;
		if (i > 0) {
			auto RanksLower = rankBox(C.rank, D.rank, i - 1);
			rank_range = RanksLower.first;
		}
	}

	///////////////////////////////////////////////////////////////////////
	///Forms the block diagonal matrix with k many copies of a, and then applies the permutations left^{-1} * block *right 

	///Significant performance improvement over the naive implementation by forming left^{-1}*block in one step
	///////////////////////////////////////////////////////////////////////////
	template<typename Derived, typename S>
	Derived permutation_block(const Eigen::PermutationMatrix<-1, -1, S>& left, const Eigen::PermutationMatrix<-1, -1, S>& right, const Eigen::MatrixBase<Derived>& a, int k) {
		int n = a.cols();
		int m = a.rows();
		Derived b = Eigen::MatrixBase<Derived>::Zero(m * k, n * k);
		for (int i = 0; i < b.rows(); i++) {
			auto u = left.indices()[i] % m;
			auto v = left.indices()[i] / m;
			b.row(i).segment(n * v, n) = a.row(u).segment(0, n);
		}
		return b * right;
	}



	/////////////////////////////////////////////////////////////////////////////
///Mixes the left and right differentials to produce the differential for the box product

///The left differentials are stored in a vector L and the right differentials in a vector R.
///The mixing is L[0], then R[0] directly to the right, then L[1] directly below then...
///It's assumed that L,R have consistent dimensions for the blocks to fit coherently.
////////////////////////////////////////////////////////////////////////////
	template<typename T>
	T MatrixMixer(std::vector<T>& L, std::vector<T>& R, int rows, int cols) {
		T mixed;
		mixed.setZero(rows, cols);
		int horz = 0, vert = 0;
		for (size_t i = 0;i < L.size();i++) {
			auto rowsL = L[i].rows();
			auto colsL = L[i].cols();
			if (rowsL > 0) {
				mixed.block(vert, horz, rowsL, colsL) = L[i];
				vert += rowsL;
			}
			auto rowsR = R[i].rows();
			auto colsR = R[i].cols();

			if (colsR > 0) {
				mixed.block(vert, horz, rowsR, colsR) = R[i];
				horz += colsR;
			}
		}
		return mixed;
	}


	///Does permutation_block and matrix mixing in one step for sparse matrices, using the triplets format.
	template<typename T, typename S, typename storage, typename U, typename V>
	void box_sparse(triplets<T, storage>& trip, const S& left, const S& right, const triplets<T, storage>& a, U rows, U cols, V copies, storage h_offset, storage v_offset) {
		for (storage k = 0; k < a.size(); k++) {
			auto i = a[k].row();
			auto j = a[k].col();
			while (i < rows * copies && j < cols * copies) {
				trip.push_back(Eigen::Triplet<T, storage>(left.indices()[i] + v_offset, right.indices()[j] + h_offset, a[k].value()));
				i += rows;
				j += cols;
			}
		}
	}


	template<typename rank_t, typename diff_t>
	void BoxPoint<rank_t, diff_t> ::getdiff() {

		if constexpr (SFINAE::is_Dense<diff_t>::value) {//dense implementation

			std::vector<diff_t> LeftDiff(i + 1);
			std::vector<diff_t> RightDiff(i + 1);
			for (int j = lowlimit; j <= highlimit; j++) {
				auto Domain = ChangeBasis<int>(C.rank[j], D.rank[i - j]);
				if (j >= 1) //We have a LeftDiff
				{
					auto RangeL = ChangeBasis<int>(C.rank[j - 1], D.rank[i - j]);
					LeftDiff[j] = permutation_block(RangeL.LefttoCanon, Domain.LefttoCanon, C.diff[j], summation(D.rank[i - j]));
				}
				if (i - j >= 1) //We have a RightDiff
				{
					typename diff_t::Scalar sign = (1 - 2 * (j % 2)); //(1 - 2 * (j % 2)) is(-1) ^ j
					auto RangeR = ChangeBasis<int>(C.rank[j], D.rank[i - j - 1]);
					RightDiff[j] = permutation_block(RangeR.RighttoCanon, Domain.RighttoCanon, static_cast<diff_t>(sign * D.diff[i - j]), summation(C.rank[j]));
				}
			}
			diff = MatrixMixer(LeftDiff, RightDiff, summation(rank_range), summation(rank_domain));

		}
		else {//sparse implementation
			typename diff_t::StorageIndex totalsize, v_offset, h_offset;
			totalsize = v_offset = h_offset = 0;

			//first find the size of nonzeros
			for (int j = lowlimit; j <= highlimit; j++) {
				if (j >= 1) 
					totalsize += C.diff[j].nonZeros()* summation(D.rank[i - j]);
				if (i - j >= 1) 
					totalsize += D.diff[i - j].nonZeros() * summation(C.rank[j]);
			}

			triplets<Scalar_t<diff_t>, typename diff_t::StorageIndex> mat;
			mat.reserve(totalsize);

			//put the nonzeros in mat
			for (int j = lowlimit; j <= highlimit; j++) {
				auto Domain = ChangeBasis<typename diff_t::StorageIndex>(C.rank[j], D.rank[i - j], 0);
				if (j >= 1) //We have a LeftDiff
				{
					auto RangeL = ChangeBasis<typename diff_t::StorageIndex>(C.rank[j - 1], D.rank[i - j], 0);
					box_sparse(mat, RangeL.LefttoCanon, Domain.LefttoCanon, make_triplets(C.diff[j]), C.diff[j].rows(), C.diff[j].cols(), summation(D.rank[i - j]), h_offset, v_offset);
					v_offset += C.diff[j].rows() * summation(D.rank[i - j]);
				}
				if (i - j >= 1) //We have a RightDiff
				{
					typename diff_t::Scalar sign = (1 - 2 * (j % 2)); //(1 - 2 * (j % 2)) is(-1) ^ j
					auto RangeR = ChangeBasis<typename diff_t::StorageIndex>(C.rank[j], D.rank[i - j - 1], 0);
					box_sparse(mat, RangeR.RighttoCanon, Domain.RighttoCanon, make_triplets(static_cast<diff_t>(sign * D.diff[i - j])), D.diff[i - j].rows(), D.diff[i - j].cols(), summation(C.rank[j]), h_offset, v_offset);
					h_offset += D.diff[i - j].cols() * summation(C.rank[j]);
				}
			}
			//turn mat to a matrix
			diff.resize(summation(rank_range), summation(rank_domain));
			diff.setFromTriplets(mat.begin(), mat.end());
		}
	}

	/// The box (tensor) product of Chains
	template<typename rank_t, typename diff_t>
	class ChainsBox :public Chains<rank_t, diff_t> {
	public:

		/// The detailed rank retains the layout of the tensor product. That's used for forming products of generators (padding)
		std::vector<std::vector<long>> detailedrank;

		///Default Constructor
		ChainsBox() {};

		/// Set the variables directly
		ChainsBox(const std::vector<rank_t>& rank, const std::vector<diff_t>& diff, const std::vector<std::vector<long>>& detailedrank)
			:Chains<rank_t, diff_t>(rank, diff), detailedrank(detailedrank) {}

		/// Get C box D up to index i.
		ChainsBox(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int);

		/// Get C box D.
		ChainsBox(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&);

	};

	template<typename rank_t, typename diff_t>
	ChainsBox<rank_t, diff_t>::ChainsBox(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int i) {
		i = std::min(i, C.maxindex + D.maxindex);
		std::vector<rank_t> rank;
		std::vector<diff_t> diff;
		std::vector<std::vector<long>> detailedrank;
		rank.reserve(i + 1);
		diff.reserve(i + 1);
		detailedrank.reserve(i + 1);
		for (int j = 0; j <= i;j++) {
			BoxPoint<rank_t, diff_t> BoxedAtPointj(C, D, j);
			rank.push_back(BoxedAtPointj.rank_domain);
			diff.push_back(BoxedAtPointj.diff);
			detailedrank.push_back(BoxedAtPointj.detailedrank);
		}
		*this = ChainsBox<rank_t, diff_t>(rank, diff, detailedrank);
	}

	template<typename rank_t, typename diff_t>
	ChainsBox<rank_t, diff_t>::ChainsBox(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D) : ChainsBox<rank_t, diff_t>(C, D, C.maxindex + D.maxindex) {}



	/// The box (tensor) product of Chains as a Junction in the given index
	template<typename rank_t, typename diff_t>
	class JunctionBox :public Junction<rank_t, diff_t> {
	public:

		/// The detailed rank retains the layout of the tensor product. That's used for forming products of generators (padding)
		std::vector<long> detailedrank;

		///Default Constructor
		JunctionBox() {};

		/// Get the i-th Junction of C box D.
		JunctionBox(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int);

		///Get the i-th Junction of an already boxed chain
		JunctionBox(const ChainsBox<rank_t, diff_t>&, int);

	};

	template<typename rank_t, typename diff_t>
	JunctionBox<rank_t, diff_t>::JunctionBox(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int i)
	{
		BoxPoint<rank_t, diff_t> OUT(C, D, i);
		detailedrank = OUT.detailedrank;
		this->rank = OUT.rank_domain;
		if (i > 0) {
			this->rankOut = std::move(OUT.rank_range);
			this->diffOut = std::move(OUT.diff);
		}
		if (i < C.maxindex + D.maxindex) {
			BoxPoint<rank_t, diff_t> IN(C, D, i + 1);
			this->rankIn = std::move(IN.rank_domain);
			this->diffIn = std::move(IN.diff);
		}
	}

	template<typename rank_t, typename diff_t>
	JunctionBox<rank_t, diff_t>::JunctionBox(const ChainsBox<rank_t, diff_t>& C, int i) : Junction<rank_t, diff_t>(C, i)
	{
		detailedrank = C.detailedrank[i];
	}


	/// Get the box product of C,D up to index i
	template<typename rank_t, typename diff_t>
	inline Chains<rank_t, diff_t> Box(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int i) {
		return ChainsBox<rank_t, diff_t>(C, D, i);
	}



	/// Get the entire box product of C,D
	template<typename rank_t, typename diff_t>
	inline Chains<rank_t, diff_t> Box(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D) {
		return Box<rank_t, diff_t>(C, D, C.maxindex + D.maxindex);
	}
}
