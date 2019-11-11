#pragma once
#include "MackeyFunctor.h"
#include "Green.h"
///@file
///@brief Contains the computations for RO(G) Chains and Homology as Mackey and Green functors. 

namespace Mackey {


	///Reindexes the homological degree k=degree[0] so that it is always 0<=k<=dimension(sphere)
	template <typename T>
	inline T Reindex(T degree) {
		auto sphere = remove0th(degree.data(), degree.size());
		degree[0] = Reindex(degree[0], sphere);
		return degree;
	}

	///Inverse of \ref Reindex "Reindex"
	template <typename T>
	inline T invReindex(T degree) {
		auto sphere = remove0th(degree.data(), degree.size());
		degree[0] = invReindex(degree[0], sphere);
		return degree;
	}

	///Returns the standard Chains of the given sphere.
	template<typename rank_t, typename diff_t, typename deg_t>
	Chains<rank_t, diff_t> StandardChains(const deg_t& sphere) {
		return StandardChains<rank_t, diff_t>(dimension(sphere), sphere);
	}


	///Returns the Chains of any given sphere up to index i
	template<typename rank_t, typename diff_t, typename deg_t>
	Chains<rank_t, diff_t> ROChains(int i, const deg_t& sphere) {
		deg_t possphere = sphere;
		deg_t negsphere = sphere;
		bool flagpos, flagneg;
		flagpos = flagneg = 0;
		for (int j = 0; j < reps; j++) {
			if (sphere[j] > 0) {
				negsphere[j] = 0;
				flagpos = 1;
			}
			else if (sphere[j] < 0) {
				possphere[j] = 0;
				flagneg = 1;
			}
		}
		if (!(flagpos && flagneg)) {
			return StandardChains<rank_t, diff_t>(i, sphere);
		}
		else {
			auto pos = StandardChains<rank_t, diff_t>(possphere);
			auto neg = StandardChains<rank_t, diff_t>(negsphere);
			return Box<rank_t, diff_t>(pos, neg, i);
		}
	}

	///Returns the Chains of any given given sphere
	template<typename rank_t, typename diff_t, typename deg_t>
	Chains<rank_t, diff_t> ROChains(const deg_t& sphere) {
		return ROChains<rank_t, diff_t>(dimension(sphere), sphere);
	}




	///The Mackey functor Homology of a Junction of Mackey functors
	template<typename rank_t, typename diff_t>
	MackeyFunctor<rank_t> HomologyLevels(const Levels<Junction<rank_t, diff_t> >& J) {
		MackeyFunctor<rank_t> M;
		M.resize(power + 1);

		std::vector<Homology<rank_t, diff_t> > H(power + 1);
		typedef Eigen::Matrix<typename diff_t::Scalar, -1, 1> gen_t;

		for (int i = 0; i < power + 1; i++) {
			H[i] = Homology<rank_t, diff_t>(J.level[i]);
			M.Groups[i] = H[i].Groups;
		}
		for (int i = 0; i < power; i++) { //for each level
			if (!H[i].isZero && !H[i + 1].isZero) {

				M.Tr[i].resize(H[i].Generators.cols());
				for (int j = 0; j < H[i].Generators.cols(); j++) { //for each generator at each level
					gen_t generator = H[i].Generators.col(j);
					M.Tr[i][j] = H[i + 1].basis(transfer(generator, J.level[i].rank, J.level[i + 1].rank));
				}

				M.Res[i].resize(H[i + 1].Generators.cols());
				for (int j = 0; j < H[i + 1].Generators.cols(); j++) { //for each generator at the higher level
					gen_t generator = H[i + 1].Generators.col(j);
					M.Res[i][j] = H[i].basis(restriction(generator, J.level[i + 1].rank, J.level[i].rank));
				}

			}
		}
		for (int i = 0; i < power; i++) { //for each level
			if (!H[i].isZero) {
				M.Weyl[i].resize(H[i].Generators.cols());
				for (int j = 0; j < H[i].Generators.cols(); j++) { //for each generator at each level
					gen_t generator = H[i].Generators.col(j);
					M.Weyl[i][j] = H[i].basis(action(generator, J.level[i].rank));
				}
			}
		}
		return M;
	}

	///Computes the Mackey functor homology of the given sphere
	template<typename rank_t, typename diff_t, typename deg_t>
	std::vector<MackeyFunctor<rank_t>> ROHomology(const deg_t& sphere)
	{
		std::vector<MackeyFunctor<rank_t>> M;
		Chains<rank_t, diff_t> A = ROChains<rank_t, diff_t>(sphere);
		for (int k = 0; k <= A.maxindex; k++) {
			Junction<rank_t, diff_t> J(A, k);
			Levels<Junction<rank_t, diff_t> > L(J);
			M.push_back(HomologyLevels<rank_t, diff_t>(L));
		}
		return M;
	}

	///Computes the Mackey functor homology of the given sphere
	template<typename rank_t, typename diff_t, typename deg_t>
	rank_t ROGreen(int level, const deg_t& first, const deg_t& second, int selectFirst, int selectSecond) {
		auto refirst = Reindex(first);
		auto resecond = Reindex(second);
		auto degreeC = refirst[0];
		int degreeD = resecond[0];
		auto firstsphere = remove0th(refirst.data(), refirst.size());
		auto secondsphere = remove0th(resecond.data(), resecond.size());
		auto limitC = std::min(degreeC + degreeD + 1, dimension(firstsphere));
		auto limitD = std::min(degreeC + degreeD + 1, dimension(secondsphere));
		Chains<rank_t, diff_t> C = ROChains<rank_t, diff_t>(limitC, firstsphere);
		Chains<rank_t, diff_t> D = ROChains<rank_t, diff_t>(limitD, secondsphere);
		Green <rank_t, diff_t>G(C, D, level, degreeC, degreeD, selectFirst, selectSecond);
		return G.normalBasis;
	}

	///Computes the Mackey functor homology of the given sphere with default 0,0 selections
	template<typename rank_t, typename diff_t, typename deg_t>
	rank_t ROGreen(int level, const deg_t& first, const deg_t& second) {
		return ROGreen<rank_t, diff_t, deg_t>(level, first, second, 0, 0);
	}

}


