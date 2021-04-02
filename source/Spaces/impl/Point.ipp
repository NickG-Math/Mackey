#pragma once
#include "../Point.hpp"

///@file
///@brief Contains the computations for RO(G) Chains and Homology as Mackey and Green functors. 

namespace mackey {

	template<typename T, typename deg_t>
	int dimension(const deg_t& sphere) {
		int sum = 0;
		for (size_t i = 0; i < T::sphere_dimensions.size(); i++)
			sum += abs(sphere[i]) * T::sphere_dimensions[i];
		return sum;
	}

	template<typename T, typename deg_t>
	deg_t Reindex(deg_t degree) {
		for (size_t i = 0; i < T::sphere_dimensions.size(); i++) {
			if (degree[i + 1] < 0)
				degree[0] -= degree[i + 1] * T::sphere_dimensions[i];
		}
		return degree;
	}

	template<typename T, typename deg_t>
	deg_t invReindex(deg_t degree) {
		for (size_t i = 0; i < T::sphere_dimensions.size(); i++) {
			if (degree[i + 1] < 0)
				degree[0] += degree[i + 1] * T::sphere_dimensions[i];
		}
		return degree;
	}

	template<typename T, typename deg_t>
	auto getChains(const deg_t& sphere, int i) {
		if (i == -1)
			i = dimension<T>(sphere);
		bool cohomology = 0;
		for (const auto& ind : sphere) {
			if (ind < 0) {
				cohomology = 1;
				break;
			}
		}
		if (!cohomology)
			return T::PositiveChains(sphere, i);
		else {
			auto C = T::PositiveChains(-sphere, dimension<T>(-sphere));
			return C.dualize(i);
		}
	}

	template<typename group_t, typename deg_t>
	auto ROChains(const deg_t& sphere, int i) {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;
		if (i == -1)
			i = dimension<group_t>(sphere);
		deg_t possphere = sphere;
		deg_t negsphere = sphere;
		bool flagpos, flagneg;
		flagpos = flagneg = 0;
		for (size_t j = 0; j < group_t::sphere_dimensions.size(); j++) {
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
			return getChains<group_t>(sphere, i);
		}
		else {
			auto pos = getChains<group_t>(possphere);
			auto neg = getChains<group_t>(negsphere);
			auto CR = Tensor<Chains<rank_t, diff_t>>(pos, neg, i).tensor();
			if constexpr (SFINAE::is_Sparse<diff_t>::value) {
				EquivariantAMT<rank_t, diff_t> Q(CR, 0, 0);
				//std::cout << Q.reduction_ratio() << "\n";
			}
			return CR;
		}
	}

	template<typename group_t>
	auto HomologyLevels(const Levels<Junction<typename group_t::rank_t, typename group_t::diff_t>, group_t >& J) {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;
		MackeyFunctor<rank_t> M;
		M.resize(group_t::power + 1);

		std::vector<Homology<rank_t, diff_t> > H(group_t::power + 1);

		for (int i = 0; i < group_t::power + 1; i++) {
			H[i] = Homology<rank_t, diff_t>(J.level[i]);
			M.Groups[i] = H[i].Groups;
		}
		for (int i = 0; i < group_t::power; i++) { //for each level
			if (!H[i].isZero && !H[i + 1].isZero) {
				M.Tr[i] = transfer<group_t>(H[i], H[i + 1], J.level[i].rank, J.level[i + 1].rank);
				M.Res[i] = restriction<group_t>(H[i + 1], H[i], J.level[i + 1].rank, J.level[i].rank);
			}
		}
		for (int i = 0; i < group_t::power; i++) { //for each level
			if (!H[i].isZero)
				M.Weyl[i] = action<group_t>(H[i], J.level[i].rank);
		}
		return M;
	}

	template<typename group_t, typename deg_t>
	auto ROHomology(int level, const deg_t& degree)
	{
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

		auto k = Reindex<group_t>(degree)[0];
		if (k < 0)
			return AbelianGroup<rank_t>(); //0 for degree reasons
		Chains<rank_t, diff_t> A = ROChains<group_t>(tail(degree.data(), degree.size(), 1));
		if (k > A.maxindex())
			return AbelianGroup<rank_t>();
		Junction<rank_t, diff_t> J(A, k);
		auto J_l = transfer<group_t>(J, level);
		Homology<rank_t, diff_t> H(J_l);
		return H.Groups;
	}

	template<typename group_t, typename deg_t>
	auto ROHomology(const deg_t& sphere)
	{
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

		std::vector<MackeyFunctor<rank_t>> M;
		Chains<rank_t, diff_t> A = ROChains<group_t>(sphere);

		for (int k = 0; k <= A.maxindex(); k++) {
			Junction<rank_t, diff_t> J(A, k);
			Levels<Junction<rank_t, diff_t>, group_t> L(J);
			M.push_back(HomologyLevels<group_t>(L));
		}
		return M;
	}

	template<typename group_t, typename deg_t>
	auto ROGreen(int level, const deg_t& first, const deg_t& second, int selectFirst, int selectSecond) {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

		auto refirst = Reindex<group_t>(first);
		auto resecond = Reindex<group_t>(second);
		auto degreeC = refirst[0];
		auto degreeD = resecond[0];
		auto firstsphere = tail(refirst.data(), refirst.size(), 1);
		auto secondsphere = tail(resecond.data(), resecond.size(), 1);
		auto limitC = std::min(degreeC + degreeD + 1, dimension<group_t>(firstsphere));
		auto limitD = std::min(degreeC + degreeD + 1, dimension<group_t>(secondsphere));
		Chains<rank_t, diff_t> C = ROChains<group_t>(firstsphere, limitC);
		Chains<rank_t, diff_t> D = ROChains<group_t>(secondsphere, limitD);
		return multiply<group_t>(C, D, level, degreeC, degreeD, selectFirst, selectSecond);
	}

	template<typename group_t, typename deg_t>
	Massey<group_t> ROMassey(int level, const deg_t& first, const deg_t& second, const deg_t& third, int selectFirst, int selectSecond, int selectThird) {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;
		auto refirst = Reindex<group_t>(first);
		auto resecond = Reindex<group_t>(second);
		auto rethird = Reindex<group_t>(third);
		auto degreeC = refirst[0];
		auto degreeD = resecond[0];
		auto degreeE = rethird[0];
		auto firstsphere = tail(refirst.data(), refirst.size(), 1);
		auto secondsphere = tail(resecond.data(), resecond.size(), 1);
		auto thirdsphere = tail(rethird.data(), rethird.size(), 1);
		Chains<rank_t, diff_t> C = ROChains<group_t>(firstsphere);
		Chains<rank_t, diff_t> D = ROChains<group_t>(secondsphere);
		Chains<rank_t, diff_t> E = ROChains<group_t>(thirdsphere);
		return Massey<group_t>(C, D, E, level, degreeC, degreeD, degreeE, selectFirst, selectSecond, selectThird, 1);
	}

}
