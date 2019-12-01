#pragma once
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <numeric>
#include <array>


///@file
///@brief Contains the construction of the change of basis matrices.
namespace Mackey {
	
	/////////////////////////////////////////////////
	/// The integer type for the permutation matrices.
	///
	/// Int for everything (including Massey products), short for everything else, char is only good for the additive structure.
	/////////////////////////////////////////////////
	typedef int pScalar; 


	///Returns the change of basis matrix (as a permutation) from basis A to basis B. A,B can be std::vector's
	template<typename Scalar, typename T>
	std::vector<Scalar> changebasis(const T& A, const T& B)
	{
		Scalar size = A.size();
		long max = *std::max_element(A.begin(), A.end());
		std::vector<Scalar> C(max + 1); //The inverse of A
		for (Scalar i = 0; i < size; i++) {
			C[A[i]] = i;
		}
		std::vector<Scalar> change;
		change.reserve(size);
		for (const auto & i : B) {
			change.push_back(C[i]); //change[i]=C[B[i]]
		}
		return change;
	}

	///Default template using the Mackey::pScalar
	template<typename Scalar=pScalar, typename T>
	std::vector<Scalar> changebasis(const T& A, const T& B);

	///Constructs the 3 bases (left+right convenient and canonical) and returns the two change of basis matrices (as permutations).
	template<typename rank_t>
	std::pair<std::vector<pScalar>, std::vector<pScalar>> boxchangebasis(const rank_t& rankC, const rank_t& rankD) {

		typedef long ScalarInternal; //Max accuracy for the ScalarInternal with 0 cost anyway

		std::vector<pScalar> LefttoCanon, RighttoCanon;

		auto s = rankC.size();
		auto t = rankD.size();
		auto totalsize = summation(rankC) * summation(rankD);

		if ((s == 1 && rankC(0) == 1) || (t == 1 && rankD(0) == 1)) { //the permutations are identities
			LefttoCanon.resize(totalsize);
			std::iota(LefttoCanon.begin(), LefttoCanon.end(), 0);
			return std::make_pair(LefttoCanon, LefttoCanon);
		}

		auto order = std::max(rankC.maxCoeff(), rankD.maxCoeff());

		std::vector<ScalarInternal> canonical;
		canonical.reserve(totalsize);

		for (int d = 0; d < t; d++) {
			for (int c = 0; c < s; c++) {
				int min, max;
				(rankC(c) < rankD(d)) ? (min = rankC(c), max = rankD(d)) : (min = rankD(d), max = rankC(c));
				for (int j = 0; j < min; j++) {
					for (int i = 0; i < max; i++) {
						canonical.push_back(c + d * s + (i % rankC(c)) * t * s + ((i + j) % rankD(d)) * order * t * s);
					}
				}
			}
		}

		std::vector<ScalarInternal> leftconv;
		leftconv.reserve(totalsize);
		for (int d = 0; d < t; d++) {
			for (int j = 0; j < rankD(d);j++) {
				for (int c = 0; c < s; c++) {
					for (int i = 0; i < rankC(c); i++) {
						leftconv.push_back(c + d * s + (i % rankC(c)) * t * s + (j % rankD(d)) * order * t * s);
					}
				}
			}
		}

		std::vector<ScalarInternal> rightconv;
		rightconv.reserve(totalsize);
		for (int c = 0; c < s; c++) {
			for (int i = 0; i < rankC(c); i++) {
				for (int d = 0; d < t; d++) {
					for (int j = 0; j < rankD(d);j++) {
						rightconv.push_back(c + d * s + (i % rankC(c)) * t * s + (j % rankD(d)) * order * t * s);
					}
				}
			}
		}
		//GCC fails with default template argument deduction in the following if I don't manually include it as <pScalar>
		return std::make_pair(changebasis<pScalar>(leftconv, canonical), changebasis<pScalar>(rightconv, canonical));
	}


	/// The change of basis permutation matrices needed for the box product
	template<typename rank_t>
	class ChangeBasis {
	public:
		Eigen::PermutationMatrix<-1, -1, pScalar> LefttoCanon; 	///<The change of basis matrix from the left convenient basis to the canonical basis
		Eigen::PermutationMatrix<-1, -1, pScalar> RighttoCanon;	///<The change of basis matrix from the right convenient basis to the canonical basis

		///Default Constructor
		ChangeBasis() {};

		/// Constructs the change of basis matrices given two ranks
		ChangeBasis(const rank_t& rank1, const rank_t& rank2) {

			typedef Eigen::Matrix<pScalar, -1, 1> vector_t;
			std::pair<std::vector<pScalar>, std::vector<pScalar>> A = boxchangebasis(rank1, rank2);

			//now just wrap into Eigen permutation matrices
			LefttoCanon = Eigen::PermutationMatrix<-1, -1, pScalar>(Eigen::Map<vector_t>(A.first.data(), A.first.size()));
			RighttoCanon = Eigen::PermutationMatrix<-1, -1, pScalar>(Eigen::Map<vector_t>(A.second.data(), A.second.size()));
		}
	};
	   	 
	////////////////////////////////////////////////////////////
/// Memoizing the ChangeBasis constructor

/// We utilize the observation that with the exception of Massey products, the ranks appearing look like ?,prime^power,...,prime^power,?. So we only need store the ?'s and total length
///////////////////////////////////////////////////////////
	template<typename rank_t>
	ChangeBasis<rank_t> memoChangeBasis(const rank_t& rank1, const rank_t& rank2) {
		static std::map<std::array<int, 6>, ChangeBasis<rank_t>> memoizer;
		const std::array<int, 6> key = {static_cast<int>(rank2.size()), rank2[0], rank2[rank2.size() - 1], rank1[0], rank1[rank1.size() - 1], static_cast<int>(rank1.size()) };
		bool notadmissible=0;
		for (int i = 1; i < rank1.size() - 1; i++) {
			if (rank1[i] != intexp(power)) {
				notadmissible = 1;
				break;
			}
		}
		for (int i = 1; i < rank2.size() - 1; i++) {
			if (rank2[i] != intexp(power)) {
				notadmissible = 1;
				break;
			}
		}
		if (!notadmissible) {
			const auto iterator = memoizer.find(key);
			if (iterator == memoizer.end()) {
				ChangeBasis<rank_t> newchange(rank1, rank2);
				memoizer[key] = newchange;
				return newchange;
			}
			else
				return iterator->second;
		}
		return ChangeBasis<rank_t>(rank1, rank2);
	}


	//A slower memoizer
	//template<typename rank_t, typename pScalar>
	//ChangeBasis<rank_t, pScalar> memoChangeBasis(const rank_t& rank1, const rank_t& rank2) {
	//	static std::map<std::pair<std::vector<typename rank_t::Scalar>, std::vector<typename rank_t::Scalar>>, ChangeBasis<rank_t>> memoizer;
	//	const std::vector<typename rank_t::pScalar> first(rank1.data(), rank1.data() + rank1.size());
	//	const std::vector<typename rank_t::pScalar> second(rank2.data(), rank2.data() + rank2.size());
	//	auto iterator = memoizer.find(std::make_pair(first, second));
	//	if (iterator == memoizer.end()) {
	//		ChangeBasis<rank_t> newchange(rank1, rank2);
	//		memoizer[std::make_pair(first, second)] = newchange;
	//		return newchange;
	//	}
	//	else {
	//		return iterator->second;
	//	}
	//}

}
