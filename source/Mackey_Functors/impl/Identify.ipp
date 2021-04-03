#pragma once
#include "../Identify.hpp"
///@file
///@brief Contains the identification classes


namespace mackey {

	///Standard equality tests if all group, group_lower, tr, res are equal 
	template<typename rank_t>
	bool IDGenerators<rank_t>::operator==(const IDGenerators<rank_t>& J) const {
		return group == J.group && group_lower == group_lower && tr == J.tr && res == J.res;
	}


	template<typename rank_t>
	MackeyFunctor<rank_t> IDGenerators<rank_t>::getMackey() const {
		MackeyFunctor<rank_t> M;
		M.resize(2);
		M.levels[0] = group_lower;
		M.levels[1] = group;
		M.tr[0] = tr;
		M.res[0] = res;
		return M;
	}


	//Use the ID data to identify all possible candidates for an element
	template<typename rank_t>
	std::vector<rank_t> id_candidates(const rank_t& basis, const IDGenerators<rank_t>& Id_first, const IDGenerators<rank_t>& Id_second) {
		auto Mack_first = Id_first.getMackey();
		auto Mack_second = Id_second.getMackey();
		auto autom = Mack_first.automorphisms();
		std::vector<size_t> match;
		match.reserve(autom.first.size());
		for (size_t i = 0; i < autom.first.size(); i++) {
			if (Mack_second == Mack_first.apply(autom.first[i], autom.second[i]))
				match.push_back(i);
		}
		std::vector<rank_t> candidates;
		candidates.reserve(match.size());
		for (auto i:match) {
			rank_t newbasis = autom.first[i][1] * basis.transpose();
			Id_second.group.normalize(newbasis);
			bool flag = 0;
			for (auto& cand : candidates) {
				if (cand == newbasis) {
					flag = 1;
					break;
				}
			}
			if (!flag)
				candidates.push_back(newbasis);
		}
		return candidates;
	}


	//Find if given elements can all be distinguished i.e. they have pairwise disjoint candidate sets
	template<typename rank_t>
	bool distinguish(const std::vector<rank_t>& elements, const IDGenerators<rank_t>& Id) {

		if (Id.group.number_of_summands() == 1) {
			for (size_t i = 0; i < elements.size(); i++) {
				for (size_t j = i + 1; j < elements.size(); j++) {
					if (elements[i] == elements[j])
						return 0;
				}
			}
			return 1;
		}

		std::vector<std::vector<rank_t>> all_candidates;
		all_candidates.reserve(elements.size());
		for (const auto& i : elements)
			all_candidates.push_back(id_candidates(i, Id, Id));
		for (size_t i = 0; i < all_candidates.size(); i++) {
			for (size_t j = i + 1; j < all_candidates.size(); j++) {
				if (std::find(all_candidates[i].begin(), all_candidates[i].end(), all_candidates[j][0]) != all_candidates[i].end())
					return 0;
			}
		}
		return 1;
	}

	namespace internal {

		//Computes the identification data given level, Junction at the bottom, and optionally if we want to store the Q matrix.
		template<typename group_t>
		IDGeneratorCompute<group_t>::IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom, bool getQ) {
			auto current = transfer<group_t>(bottom, level);
			rank_level = current.rank;
			H_level = Homology<rank_t, diff_t>(current, getQ);

			ID.group = H_level.group;
			if (ID.group.iscyclic())
				return;
			else
				id(transfer<group_t>(bottom, level - 1));
		}

		//Identifies the generators in case of non cyclic homology.
		template<typename group_t>
		void IDGeneratorCompute<group_t>::id(const Junction<rank_t, diff_t>& lower) {
			auto rank_lower = lower.rank;
			Homology<rank_t, diff_t>H_lower(lower);
			if (!H_lower.isZero) {
				ID.tr = transfer<group_t>(H_lower, H_level, rank_lower, rank_level);
				ID.res = restriction<group_t>(H_level, H_lower, rank_level, rank_lower);
				ID.group_lower = H_lower.group;
			}
		}
	}

}
