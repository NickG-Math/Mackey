#pragma once
#include "../Mult_Identify.hpp"

///@file
///@brief Contains extra identification methods for the multiplication graph.

namespace mackey {

	template<typename group_t>
	MultIdentify<group_t>::MultIdentify(MultGraph<group_t>& MG) : MG(MG) {
		triples_to_be_done.reserve(MG.unidentified.size());
	}

	template<typename group_t>
	void MultIdentify<group_t>::pass_all_unidentified() {
		int size = 0;
		while (size != MG.unidentified.size() && MG.unidentified.size() > 0) {
			size = MG.unidentified.size();
			for (const auto& pair : MG.unidentified)
				identify_connection(pair.first, pair.second); //this adds extra unidentified's but no more triples so it's fine to stop at size.
			MACKEY_RUN_LOOP_PARALLEL
				for (int i = 0; i < triples_to_be_done.size(); i++) {
					auto G = MG.triple_product(triples_to_be_done[i][0], triples_to_be_done[i][1], triples_to_be_done[i][2]);
					MACKEY_RUN_BLOCK_SERIAL
						MG.tripleGreens[{triples_to_be_done[i][0], triples_to_be_done[i][1], triples_to_be_done[i][2]}] = G;
				}
			triples_to_be_done.clear();
			for (const auto& triple : identified) {
				make_connection(triple[0], triple[1], triple[2]);
				MG.unidentified.erase(std::pair(triple[0],triple[1]));
			}
			identified.clear();
		}
	}

	template<typename group_t>
	void MultIdentify<group_t>::identify_connection(int i, int j1) {
		auto prod = MG.product_element_irreducible(i, j1);
		auto candidates = MG.product_candidates(prod);
		for (const auto& cand : candidates) { //add all candidates and their edges
			int t = MG.find_and_add_element(prod.deg, cand);
			if (t + 1 != MG.element.size()) //nothing added
				continue;
			MG.add_edges(t);
		}
		auto distinguished = distinguish(prod.deg, candidates); //see if candidates can be distinguished from the edges out of them
		int j2 = distinguished.first;
		if (j2 == -1)
			return;
		identified.insert({ i,j1,j2 });
		int deg_i = MG.tracker[i];
		if (MG.tripleGreens.count({ deg_i,j1,j2 }) == 0 && std::find(triples_to_be_done.begin(), triples_to_be_done.end(), std::array{ deg_i,j1,j2 }) == triples_to_be_done.end())
			triples_to_be_done.push_back({ deg_i,j1,j2 });
	}


	template<typename group_t>
	void MultIdentify<group_t>::make_connection(int i, int j1, int j2) {
		auto prod = MG.product_element_irreducible(i, j1);
		auto candidates = MG.product_candidates(prod);
		auto distinguished = distinguish(prod.deg, candidates, j2);
		auto prod_all = MG.product_element_irreducible(i, j1, j2);
		if (prod_all.basis.isZero() || prod_all.basis.size() == 1) {
			int k = std::find(distinguished.begin(), distinguished.end(), prod_all.basis) - distinguished.begin();
			MG.connect(candidates[k], i, j1, prod.deg);
		}
		else {
			auto cand_prod = MG.product_candidates(prod_all);
			for (const auto& k : distinguished) {
				auto it = std::find(cand_prod.begin(), cand_prod.end(), k);
				if (it != cand_prod.end())
					MG.connect(candidates[it - cand_prod.begin()], i, j1, prod.deg);
			}
		}
	}


	template<typename group_t>
	auto MultIdentify<group_t>::distinguish(int degree_ind, const std::vector<rank_t>& candidates) -> std::pair<int, std::vector<rank_t>> {

		for (int j = 0; j < MG.basicIrreducibles.size(); j++) { //for each basic irreducible
			auto prods = distinguish(degree_ind, candidates, j);
			if (!prods.empty())
				return std::pair(j, prods);
		}
		return std::pair(-1, std::vector<rank_t>());
	}

	template<typename group_t>
	auto MultIdentify<group_t>::distinguish(int degree_ind, const std::vector<rank_t>& candidates, int j) -> std::vector<rank_t> {

		const auto degree_prod = MG.index_product(degree_ind, j);
		if (degree_prod == -1)
			return std::vector<rank_t>();

		const auto Id = MG.NonZeroHomology.at(MG.degree[degree_prod]);
		const rank_t zero_matrix = rank_t::Zero(Id.group.number_of_summands());

		std::vector<rank_t> products;
		products.reserve(candidates.size());

		bool inadmissible = 0;

		for (const auto& cand : candidates) { //for each candidate
			int i = MG.antielement.at(std::make_pair(degree_ind, cand)); //get element index

			bool found = 0;
			for (auto it = MG.graph.neighbors(i).begin(); it != MG.graph.neighbors(i).end(); ++it) //search to see if i*j is known
				if (it.edge().id == j && it.edge().color == 0) { // success!
					products.push_back(MG.element[it.node()]); // products.push_back(MG.product_element_irreducible(i, j).basis);
					found = 1;
					break;
				}
			if (!found) { //fail: i multiplies with j to something unknown or to 0
				if (std::find(MG.zeroproduct[i].begin(), MG.zeroproduct[i].end(), j) == MG.zeroproduct[i].end()) {
					inadmissible = 1; //it's not 0
					break;
				}
				else
					products.push_back(zero_matrix); //product is 0
			}
		}
		if (inadmissible || !mackey::distinguish(products, Id)) //if you can't distinguish
			return std::vector<rank_t>();
		else //if you can distinguish
			return products;
	}

}

