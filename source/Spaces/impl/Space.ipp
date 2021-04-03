#pragma once
#include "../Space.hpp"

///@file
///@brief File containing the class for spaces other than points!

namespace mackey {

	template<typename space_t, typename group_t>
	const auto& Space<space_t, group_t>::getChains() {
		if (complex.diff.empty())
			setcomplexes();
		return complex;
	}

	template<typename space_t, typename group_t>
	const auto& Space<space_t, group_t>::getCoChains() {
		if (cocomplex.diff.empty())
			setcomplexes();
		return cocomplex;
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::ROHomology(int level, const std::vector<int>& homologysphere) {
		auto C = getROChains(homologysphere);
		std::vector<typename group_t::rank_t> groups;
		groups.reserve(C.maxindex + 1);
		for (int k = 0; k <= C.maxindex; k++) {
			Junction<rank_t, diff_t> J(C, k);
			auto J_l = transfer<group_t>(J, level);
			Homology<rank_t, diff_t> H(J_l);
			groups.push_back(H.group);
		}
		return groups;
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::ROHomology(const std::vector<int>& homologysphere) {
		auto C = getROChains(homologysphere);
		std::vector<MackeyFunctor<typename group_t::rank_t>> M;
		M.reserve(C.maxindex + 1);
		for (int k = 0; k <= C.maxindex; k++) {
			Junction<rank_t, diff_t> J(C, k);
			Levels<Junction<rank_t, diff_t>, group_t> L(J);
			auto S = HomologyLevels<group_t>(L);
			S.notation();
			M.push_back(S);
		}
		return M;
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::ROCohomology(int level, const std::vector<int>& cohomologysphere) {
		auto C = getROCoChains(cohomologysphere);
		std::vector<rank_t> groups;
		groups.reserve(C.maxindex + 1);
		for (int k = C.maxindex; k >= 0; k--) {
			Junction<rank_t, diff_t> J(C, k);
			auto J_l = transfer<group_t>(J, level);
			Homology<rank_t, diff_t> H(J_l);
			groups.push_back(H.group);
		}
		return groups;
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::ROCohomology(const std::vector<int>& cohomologysphere) {
		auto C = getROCoChains(cohomologysphere);
		std::vector<MackeyFunctor<rank_t>> M;
		M.reserve(C.maxindex() + 1);
		for (int k = C.maxindex(); k >= 0; k--) {
			Junction<rank_t, diff_t> J(C, k);
			Levels<Junction<rank_t, diff_t>, group_t> L(J);
			auto S = HomologyLevels<group_t>(L);
			S.notation();
			M.push_back(S);
		}
		return M;
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::boundary(int dimension_domain, int cell) {
		return static_cast<space_t*>(this)->boundary(dimension_domain, cell);
	}

	template<typename space_t, typename group_t>
	long Space<space_t, group_t>::estimate_nonzero_entries(int dimension_domain) {
		return static_cast<space_t*>(this)->estimate_nonzero_entries(dimension_domain);
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::getROChains(const std::vector<int>& sphere) {
		auto C = ROChains<group_t>(-sphere);
		return Tensor<Chains<rank_t, diff_t>>(C, complex).tensor();
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::getROCoChains(const std::vector<int>& sphere) {
		auto C = ROChains<group_t>(-sphere).dualize();
		return Tensor<Chains<rank_t, diff_t>>(C, cocomplex).tensor();
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::getdiff(int domain) {
		const rank_t& rank_domain = cells[domain];
		const rank_t& rank_range = cells[domain - 1];
		diff_t diff;
		std::conditional_t<SFINAE::is_Sparse<diff_t>::value, std::vector<triplet_t<diff_t>>, int> triplets; //only used when sparse
		if constexpr (SFINAE::is_Dense<diff_t>::value)
			diff = diff_t::Zero(summation(rank_range), summation(rank_domain));
		else {
			diff = diff_t(summation(rank_range), summation(rank_domain));
			auto i = estimate_nonzero_entries(domain);
			if (i >= 0)
				triplets.reserve(i);
			else
				triplets.reserve(summation(rank_range));
		}
		for (int j = 0; j < diff.cols(); j++) {
			auto bd = boundary(domain, j);
			for (auto& o : bd) {
				if constexpr (SFINAE::is_Dense<diff_t>::value)
					diff(o.cell, j) = o.coefficient;
				else
					triplets.emplace_back(o.cell, j, o.coefficient);
			}
		}
		if constexpr (SFINAE::is_Sparse<diff_t>::value)
			diff.setFromTriplets(triplets.begin(), triplets.end());
		return diff;
	}

	template<typename space_t, typename group_t>
	void Space<space_t, group_t>::setcomplexes() {
		std::vector<diff_t> diff;
		diff.reserve(cells.size());
		diff.push_back(diff_t());
		for (int i = 0; i < cells.size() - 1; i++)
			diff.push_back(getdiff(i + 1));
		complex = Chains<rank_t, diff_t>(cells, diff);
		cocomplex = complex.dualize();
	}

}
