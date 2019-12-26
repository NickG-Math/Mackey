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
	/// Use through JunctionBox or the friend function Box.
	//////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	class BoxPoint {
		rank_t rank_domain, rank_range;
		diff_t diff;

		/// The detailed rank retains the layout of the tensor product. That's used for forming products of generators (padding)
		std::vector<rank_t> detailedrank;

		BoxPoint(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int i)
			: C(C), D(D), i(i), lowlimit(std::max(0, i - D.maxindex)), highlimit(std::min(i, C.maxindex)) {}
		const Chains<rank_t, diff_t>& C, D;
		const int i, lowlimit, highlimit;
		void getranks();
		void getdomainrank();
		void getrangerank();
		void getdiff();

		/// The box (tensor) product of Chains as a Junction in a desired degree
		template<typename shadow_rank_t, typename shadow_diff_t>
		friend class JunctionBox;

		/// The box (tensor) product of Chains
		template<typename shadow_rank_t, typename shadow_diff_t>
		friend class ChainsBox;

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
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < B.size(); j++) {
				C.segment(current, std::min(A(i), B(j))).setConstant(std::max(A(i), B(j)));
				current += std::min(A(i), B(j));
			}
		}
		return C;
	}

	///The rank of the i-th differential of the tensor product of chains C,D
	template<typename rank_t, typename diff_t>
	std::pair<rank_t, std::vector<rank_t>> rankBox(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int i) {
		rank_t rank;
		std::vector<rank_t> detailedrank;
		detailedrank.resize(i + 1);
		auto lowlimit = std::max(0, i - D.maxindex);
		auto highlimit = std::min(i, C.maxindex);
		int sum = 0;
		for (int j = lowlimit; j <= highlimit; j++) {
			detailedrank[j] = rankmult(C.rank[j], D.rank[i - j]);
			sum += detailedrank[j].size();
		}
		rank.resize(sum);
		int track = 0;
		for (const auto& i : detailedrank) {
			rank.segment(track, i.size()) = i;
			track += i.size();
		}
		return std::make_pair(rank, detailedrank);
	}

	template<typename rank_t, typename diff_t>
	void BoxPoint<rank_t, diff_t> ::getranks() {
		getdomainrank();
		if (i > 0) {
			getrangerank();
		}
	}

	template<typename rank_t, typename diff_t>
	void BoxPoint<rank_t, diff_t> ::getdomainrank() {
		auto Ranks = rankBox(C, D, i);
		rank_domain = Ranks.first;
		detailedrank = Ranks.second;
	}

	template<typename rank_t, typename diff_t>
	void BoxPoint<rank_t, diff_t> ::getrangerank() {
		auto RanksLower = rankBox(C, D, i - 1);
		rank_range = RanksLower.first;
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
	template<typename Derived>
	Derived MatrixMixer(std::vector<Derived>& L, std::vector<Derived>& R) {
		auto rows = L[0].rows();
		auto cols = std::max(L[0].cols(), R[0].cols());
		for (size_t i = 0; i < L.size() - 1; i++) {
			rows += std::max(L[i + 1].rows(), R[i].rows());
			cols += std::max(L[i + 1].cols(), R[i + 1].cols());
		}
		rows += R[R.size() - 1].rows();

		Derived mixed;
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



	template<typename rank_t, typename diff_t>
	void BoxPoint<rank_t, diff_t> ::getdiff() {
		std::vector<diff_t> LeftDiff(i + 1);
		std::vector<diff_t> RightDiff(i + 1);

		for (int j = lowlimit; j <= highlimit; j++) {
			auto Domain = memoChangeBasis<rank_t>(C.rank[j], D.rank[i - j]);
			if (j >= 1) //We have a LeftDiff
			{
				auto RangeL = memoChangeBasis<rank_t>(C.rank[j - 1], D.rank[i - j]);
				LeftDiff[j] = permutation_block(RangeL.LefttoCanon, Domain.LefttoCanon, C.diff[j], summation(D.rank[i - j]));
			}
			if (i - j >= 1) //We have a RightDiff
			{
				typename diff_t::Scalar sign = (1 - 2 * (j % 2)); //(1 - 2 * (j % 2)) is(-1) ^ j
				auto RangeR = memoChangeBasis<rank_t>(C.rank[j], D.rank[i - j - 1]);
				RightDiff[j] = permutation_block(RangeR.RighttoCanon, Domain.RighttoCanon, static_cast<diff_t>(sign * D.diff[i - j]), summation(C.rank[j]));
			}
		}
		diff = MatrixMixer(LeftDiff, RightDiff);
	}

	/// The box (tensor) product of Chains
	template<typename rank_t, typename diff_t>
	class ChainsBox :public Chains<rank_t, diff_t> {
	public:

		/// The detailed rank retains the layout of the tensor product. That's used for forming products of generators (padding)
		std::vector<std::vector<rank_t>> detailedrank;

		///Default Constructor
		ChainsBox() {};

		/// Set the variables directly
		ChainsBox(const std::vector<rank_t>& rank, const std::vector<diff_t>& diff, const std::vector<std::vector<rank_t>>& detailedrank)
			:Chains<rank_t, diff_t>(rank, diff), detailedrank(detailedrank) {}

		/// Get C box D up to index i.
		ChainsBox(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int);

		/// Get C box D.
		ChainsBox(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&);

	};

	template<typename rank_t, typename diff_t>
	ChainsBox<rank_t, diff_t>::ChainsBox(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int i) {
		std::vector<rank_t> rank;
		std::vector<diff_t> diff;
		std::vector<std::vector<rank_t>> detailedrank;
		rank.reserve(i + 1);
		diff.reserve(i + 1);
		detailedrank.reserve(i + 1);
		for (int j = 0; j <= i;j++) {
			BoxPoint<rank_t, diff_t> BoxedAtPointj(C, D, j);
			BoxedAtPointj.getdomainrank();
			BoxedAtPointj.getdiff();
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
		std::vector<rank_t> detailedrank;

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
		OUT.getranks();
		detailedrank = OUT.detailedrank;
		this->rank = OUT.rank_domain;
		if (i > 0) {
			OUT.getdiff();
			this->rankOut = std::move(OUT.rank_range);
			this->diffOut = std::move(OUT.diff);
		}
		if (i < C.maxindex + D.maxindex) {
			BoxPoint<rank_t, diff_t> IN(C, D, i + 1);
			IN.getdomainrank();
			IN.getdiff();
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
