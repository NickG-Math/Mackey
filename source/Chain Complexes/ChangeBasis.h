#pragma once
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <numeric>
#include <array>


///@file
///@brief Contains the construction of the change of basis matrices.

namespace {
	template<typename, typename, typename = void>
	struct inv_type;

	template<typename pScalar, typename T>
	struct inv_type<pScalar, T, typename std::enable_if_t<std::is_integral<T>::value>> {
		typedef std::vector<pScalar> t;
	};

	template<typename pScalar, typename T>
	struct inv_type<pScalar, T, typename std::enable_if_t<!std::is_integral<T>::value>> {
		typedef std::map<T, pScalar> t;
	};
}

namespace Mackey {

	///Returns the change of basis matrix (as a permutation) from basis A to basis B.
	template<typename pScalar, typename T>
	std::vector<pScalar> changebasis(const std::vector<T>& A, const std::vector<T>& B)
	{
		typename inv_type<pScalar,T>::t A_inverse; //The inverse of A
		if constexpr (std::is_integral<T>::value) {
			long max = *std::max_element(A.begin(), A.end());
			A_inverse.resize(max + 1); 
		}
		for (pScalar i = 0; i <A.size(); i++)
			A_inverse[A[i]] = i;
		std::vector<pScalar> change;
		change.reserve(A.size());
		for (const auto& i : B)
			change.push_back(A_inverse[i]); //change=A_inverse B
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
			return std::make_pair(changebasis<pScalar>(leftconv, canonical), changebasis<pScalar>(rightconv, canonical));
		return std::make_pair(changebasis<pScalar>(canonical, leftconv), changebasis<pScalar>(canonical, rightconv));

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
}
