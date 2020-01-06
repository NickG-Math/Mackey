#pragma once
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <numeric>
#include <array>


///@file
///@brief Contains the construction of the change of basis matrices.
namespace Mackey {

	///Returns the change of basis matrix (as a permutation) from basis A to basis B. 
	template<typename pScalar, typename T>
	std::vector<pScalar> changebasis(const T& A, const T& B)
	{
		pScalar size = A.size();
		long max = *std::max_element(A.begin(), A.end());
		std::vector<pScalar> C(max + 1); //The inverse of A
		for (pScalar i = 0; i < size; i++)
			C[A[i]] = i;
		std::vector<pScalar> change;
		change.reserve(size);
		for (const auto& i : B)
			change.push_back(C[i]); //change[i]=C[B[i]]
		return change;
	}

	///Constructs the 3 bases (left+right convenient and canonical) and returns the two change of basis matrices (as permutations).
	template<typename pScalar, typename rank_t>
	std::pair<std::vector<pScalar>, std::vector<pScalar>> boxchangebasis(const rank_t& rankC, const rank_t& rankD, bool convtocanon) {

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
		if (convtocanon)
			return std::make_pair(changebasis<pScalar, std::vector<ScalarInternal>>(leftconv, canonical), changebasis<pScalar, std::vector<ScalarInternal>>(rightconv, canonical));
		return std::make_pair(changebasis<pScalar, std::vector<ScalarInternal>>(canonical, leftconv), changebasis<pScalar, std::vector<ScalarInternal>>(canonical, rightconv));

	}


	/// The change of basis permutation matrices needed for the box product
	template<typename pScalar>
	class ChangeBasis {
	public:
		Eigen::PermutationMatrix<-1, -1, pScalar> LefttoCanon; 	///<The change of basis matrix from the left convenient basis to the canonical basis
		Eigen::PermutationMatrix<-1, -1, pScalar> RighttoCanon;	///<The change of basis matrix from the right convenient basis to the canonical basis

		///Default constructor
		ChangeBasis() {};

		/// Constructs the change of basis matrices given two ranks. The boolean variable determines whether we go from the convenient to canonical or the opposite
		template<typename rank_t>
		ChangeBasis(const rank_t& rank1, const rank_t& rank2, bool convtocanon) {
			typedef Eigen::Matrix<pScalar, -1, 1> vector_t;
			auto A = boxchangebasis<pScalar,rank_t>(rank1, rank2, convtocanon);
			//now just wrap into Eigen permutation matrices
			LefttoCanon = Eigen::PermutationMatrix<-1, -1, pScalar>(Eigen::Map<vector_t>(A.first.data(), A.first.size()));
			RighttoCanon = Eigen::PermutationMatrix<-1, -1, pScalar>(Eigen::Map<vector_t>(A.second.data(), A.second.size()));
		}

		/// Constructs the change of basis matrices given two ranks, from the convenient bases to the canonical one.
		template<typename rank_t>
		ChangeBasis(const rank_t& rank1, const rank_t& rank2) : ChangeBasis(rank1, rank2, 1) {}

	};

//	////////////////////////////////////////////////////////////
///// Memoizing the ChangeBasis constructor. Not thread safe hence why we lock!
//
///// We utilize the observation that the ranks appearing usually look like ?,prime^power,...,prime^power,?. So we only need store the ?'s and total length
/////////////////////////////////////////////////////////////
//	template<typename rank_t>
//	ChangeBasis<int> memoChangeBasis(const rank_t& rank1, const rank_t& rank2, bool convtocanon) {
//		static std::map<std::array<int, 7>, ChangeBasis<int>> memoizer;
//		const std::array<int, 7> key = { static_cast<int>(rank2.size()), rank2[0], rank2[rank2.size() - 1], rank1[0], rank1[rank1.size() - 1], static_cast<int>(rank1.size()), convtocanon };
//		bool notadmissible = 0;
//		for (int i = 1; i < rank1.size() - 1; i++) {
//			if (rank1[i] != intexp(power)) {
//				notadmissible = 1;
//				break;
//			}
//		}
//		for (int i = 1; i < rank2.size() - 1; i++) {
//			if (rank2[i] != intexp(power)) {
//				notadmissible = 1;
//				break;
//			}
//		}
//		if (notadmissible)
//			return ChangeBasis<int>(rank1, rank2, convtocanon);
//		ChangeBasis<int> change;
//
//#ifdef PARALLELIZE
//		omp_set_lock(&Mackeylocker);
//#endif
//		const auto iterator = memoizer.find(key);
//		if (iterator == memoizer.end()) {
//			change=ChangeBasis<int>(rank1, rank2, convtocanon);
//			memoizer[key] = change;
//		}
//		else
//			change=iterator->second;
//#ifdef PARALLELIZE
//		omp_unset_lock(&Mackeylocker);
//#endif
//		return change;
//	}
//
//
//	////////////////////////////////////////////////////////////
///// Memoizing the ChangeBasis constructor. Not thread safe hence why we lock!
//
///// We utilize the observation that the ranks appearing usually look like ?,prime^power,...,prime^power,?. So we only need store the ?'s and total length
/////////////////////////////////////////////////////////////
//	template<typename rank_t>
//	ChangeBasis<int> memoChangeBasis(const rank_t& rank1, const rank_t& rank2) {
//		return memoChangeBasis(rank1, rank2, 1);
//	}

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
