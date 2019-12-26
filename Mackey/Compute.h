#pragma once
#include "MackeyFunctor.h"
#include "Green.h"
#include "Massey.h"
///@file
///@brief Contains the computations for RO(G) Chains and Homology as Mackey and Green functors. 

namespace Mackey {


	///Reindexes the homological degree k=degree[0] so that it is always 0<=k<=dimension(sphere)
	template <typename T>
	inline T Reindex(T degree) {
		auto sphere = tail(degree.data(), degree.size());
		degree[0] = Reindex(degree[0], sphere);
		return degree;
	}

	///Inverse of \ref Reindex "Reindex"
	template <typename T>
	inline T invReindex(T degree) {
		auto sphere = tail(degree.data(), degree.size());
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

		for (int i = 0; i < power + 1; i++) {
			H[i] = Homology<rank_t, diff_t>(J.level[i]);
			M.Groups[i] = H[i].Groups;
		}
		for (int i = 0; i < power; i++) { //for each level
			if (!H[i].isZero && !H[i + 1].isZero) {
				M.Tr[i] = transfer(H[i], H[i + 1], J.level[i].rank, J.level[i + 1].rank);
				M.Res[i] = restriction(H[i + 1], H[i], J.level[i + 1].rank, J.level[i].rank);
			}
		}
		for (int i = 0; i < power; i++) { //for each level
			if (!H[i].isZero)
				M.Weyl[i] = action(H[i], J.level[i].rank);
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

	///Computes the product of two generators in the RO(G) homology given their level, degrees and selections (if noncyclic)
	template<typename rank_t, typename diff_t, typename deg_t>
	rank_t ROGreen(int level, const deg_t& first, const deg_t& second, int selectFirst, int selectSecond) {
		auto refirst = Reindex(first);
		auto resecond = Reindex(second);
		auto degreeC = refirst[0];
		auto degreeD = resecond[0];
		auto firstsphere = tail(refirst.data(), refirst.size());
		auto secondsphere = tail(resecond.data(), resecond.size());
		auto limitC = std::min(degreeC + degreeD + 1, dimension(firstsphere));
		auto limitD = std::min(degreeC + degreeD + 1, dimension(secondsphere));
		Chains<rank_t, diff_t> C = ROChains<rank_t, diff_t>(limitC, firstsphere);
		Chains<rank_t, diff_t> D = ROChains<rank_t, diff_t>(limitD, secondsphere);
		return multiply(C, D, level, degreeC, degreeD, selectFirst, selectSecond);
	}

	///Computes the product of two generators in the RO(G) homology given their level, degrees and default 0,0 selections (if noncyclic)
	template<typename rank_t, typename diff_t, typename deg_t>
	inline rank_t ROGreen(int level, const deg_t& first, const deg_t& second) {
		return ROGreen<rank_t, diff_t, deg_t>(level, first, second, 0, 0);
	}

	///Computes the Massey product of three generators in the RO(G) homology given their level, degrees and selections (if noncyclic)
	template<typename rank_t, typename diff_t, typename deg_t>
	Massey<rank_t, diff_t> ROMassey(int level, const deg_t& first, const deg_t& second, const deg_t& third, int selectFirst, int selectSecond, int selectThird) {
		auto refirst = Reindex(first);
		auto resecond = Reindex(second);
		auto rethird = Reindex(third);
		auto degreeC = refirst[0];
		auto degreeD = resecond[0];
		auto degreeE = rethird[0];
		auto firstsphere = tail(refirst.data(), refirst.size());
		auto secondsphere = tail(resecond.data(), resecond.size());
		auto thirdsphere = tail(rethird.data(), rethird.size());
		Chains<rank_t, diff_t> C = ROChains<rank_t, diff_t>(firstsphere);
		Chains<rank_t, diff_t> D = ROChains<rank_t, diff_t>(secondsphere);
		Chains<rank_t, diff_t> E = ROChains<rank_t, diff_t>(thirdsphere);
		return Massey<rank_t, diff_t>(C, D, E, level, degreeC, degreeD, degreeE, selectFirst, selectSecond, selectThird);
	}

	///Computes the Massey product of three generators in the RO(G) homology given their level, degrees and default 0,0,0 selections (if noncyclic)
	template<typename rank_t, typename diff_t, typename deg_t>
	inline Massey<rank_t, diff_t> ROMassey(int level, const deg_t& first, const deg_t& second, const deg_t& third) {
		return ROMassey<rank_t, diff_t, deg_t>(level, first, second, third, 0, 0, 0);
	}


}


