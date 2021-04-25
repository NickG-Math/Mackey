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
		std::vector<AbelianGroup<rank_t>> groups;
		groups.reserve(C.maxindex() + 1);
		for (int k = 0; k <= C.maxindex(); k++) {
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
		std::vector<MackeyFunctor<rank_t>> M;
		M.reserve(C.maxindex() + 1);
		for (int k = 0; k <= C.maxindex(); k++) {
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
		std::vector<AbelianGroup<rank_t>> groups;
		groups.reserve(C.maxindex() + 1);
		for (int k = C.maxindex(); k >= 0; k--) {
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
	Space<space_t, group_t>::equivariant_index::equivariant_index(int index) :index(index) {}

	template<typename space_t, typename group_t>
	Space<space_t, group_t>::nonequivariant_index::nonequivariant_index(int index) : index(index) {}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::convert(nonequivariant_index n, int dimension) const -> equivariant_index {
		auto index = 0;
		auto sum = 0;
		for (;sum < n.index && index < cells[dimension].size(); index++)
			sum += cells[dimension][index];
		if (sum!= n.index) {
			std::cerr << "Nonequivariant index provided must come from equivariant index";
			abort();
		}
		return equivariant_index(index);
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::convert(equivariant_index e, int dimension) const -> nonequivariant_index {
		auto index = 0;
		bool flag = 0;
		for (int i = 0; i < e.index; i++) {
			index += cells[dimension][i];
		}
		return nonequivariant_index(index);
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::action(nonequivariant_index n, int dimension) const -> nonequivariant_index {
		auto start = 0;
		int i = 0;
		int length;
		for (;start < n.index && i < cells[dimension].size(); i++)
			start += cells[dimension][i];
		if (start > n.index) {
			start -= cells[dimension][i - 1];
			length = cells[dimension][i-1];
		}
		else
			length = cells[dimension][i];
		auto nextindex=n.index+1;
		if (nextindex >= start + length)
			nextindex = start;
		return nonequivariant_index(nextindex);
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::boundary(int dimension_domain, equivariant_index cell) {
		return static_cast<space_t*>(this)->boundary(dimension_domain, cell);
	}

	template<typename space_t, typename group_t>
	int64_t Space<space_t, group_t>::estimate_nonzero_entries(int dimension_domain) {
		return static_cast<space_t*>(this)->estimate_nonzero_entries(dimension_domain);
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::getROChains(const std::vector<int>& sphere) {
		setcomplexes();
		auto C = ROChains<group_t>(-sphere);
		return Tensor<Chains<rank_t, diff_t>>(C, complex).tensor();
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::getROCoChains(const std::vector<int>& sphere) {
		setcomplexes();
		auto C = ROChains<group_t>(-sphere).dualize();
		return Tensor<Chains<rank_t, diff_t>>(C, cocomplex).tensor();
	}

	template<typename space_t, typename group_t>
	auto Space<space_t, group_t>::getdiff(int domain) -> diff_t {
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
		for (int j = 0; j < rank_domain.size(); j++) {
			auto eq = equivariant_index(j);
			nonequivariant_index neq = convert(eq,domain);
			auto bd = boundary(domain, eq);
			for (auto& o : bd) {
				nonequivariant_index next_domain = neq;
				nonequivariant_index next_range= o.cell;
				for (int g = 0; g < cells[domain][j]; g++) {
					if constexpr (SFINAE::is_Dense<diff_t>::value)
						diff(next_range.index, next_domain.index) = o.coefficient;
					else
						triplets.emplace_back(next_range.index, next_domain.index, o.coefficient);
					next_domain=action(next_domain, domain);
					next_range=action(next_range, domain-1);
				}
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
