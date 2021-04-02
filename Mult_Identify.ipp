#pragma once
#include "../Mult_Identify.hpp"

///@file
///@brief Contains extra identification methods for the multiplication graph.

namespace mackey {

	template<typename group_t>
	MultIdentify<group_t>::MultIdentify(MultGraph<group_t>& MG) : MG(MG) {
		triples_to_be_done.reserve(MG.unidentified.size());
		identified.reserve(MG.unidentified.size());
		can_do_more = 1;
	}

	template<typename group_t>
	void MultIdentify<group_t>::pass_all_unidentified() {
		int size = 0;
		while (size != MG.unidentified.size() && MG.unidentified.size() > 0) {
			size = MG.unidentified.size();
			for (int i = 0; i < size; i++) {
				auto thepair = MG.unidentified[i];
				pass_triple(thepair, 1); //this adds extra unidentified's but no more triples so it's fine to stop at size.
			}
			pass_identified();
		}
	}

	template<typename group_t>
	void MultIdentify<group_t>::pass_identified() {

		MACKEY_RUN_LOOP_PARALLEL
			for (int i = 0; i < triples_to_be_done.size(); i++) {
				auto G = MG.triple_product(triples_to_be_done[i][0], triples_to_be_done[i][1], triples_to_be_done[i][2]);
				MACKEY_RUN_BLOCK_SERIAL{
				MG.tripleGreens[{triples_to_be_done[i][0], triples_to_be_done[i][1], triples_to_be_done[i][2]}] = G;
				}
			}
		triples_to_be_done.clear();
		for (const auto& pair : identified)
			pass_triple(pair, 0);
		delete_identified();
	}

	template<typename group_t>
	void MultIdentify<group_t>::delete_identified() {
		std::vector<std::pair<int, int>> true_unidentified;
		true_unidentified.reserve(MG.unidentified.size() - identified.size());
		for (auto i : MG.unidentified)
			if (std::find(identified.begin(), identified.end(), i) == identified.end())
				true_unidentified.push_back(i);
		MG.unidentified = true_unidentified;
		identified.clear();
	}

	template<typename group_t>
	bool MultIdentify<group_t>::pass_triple(std::pair<int, int> ij, bool collection_stage) {
		auto i = ij.first;
		auto j = ij.second;
		auto deg_i = MG.tracker[i];
		auto G = MG.Greens[deg_i][j];
		int deg_ij = MG.index_product(deg_i, j); //can't be -1
		auto v = determine_connection_identify(G, i, j, deg_ij, collection_stage);

		if (collection_stage && v.size() != 0)
			return 1;
		else if (v.size() != 0) {
			MG.connect(v, i, j, deg_ij);
			return 1;
		}
		return 0;
	}

	template<typename group_t>
	auto MultIdentify<group_t>::determine_connection_identify(const Green<group_t>& G_ij1, int i, int j1, int deg_ij1, bool collection_stage) -> rank_t {
		rank_t basis = G_ij1.getNormalBasis(0, MG.element[i]);

		std::vector<rank_t> candidates = id_candidates(basis, G_ij1.boxID, MG.NonZeroHomology.at(MG.degree[deg_ij1]));

		for (const auto& cand : candidates) { //add all candidates and their edges
			int t = MG.find_and_add_element(deg_ij1, cand);
			if (t+1 != MG.element.size()) //nothing added
				continue;
			MG.add_edges(t);
		}

		std::pair<int, std::vector<rank_t>> distinguished = distinguish(deg_ij1, candidates); //see if candidates can be distinguished from the edges out of them
		int j2 = distinguished.first;
		if (j2 == -1)
			return rank_t();

		if (std::find(identified.begin(), identified.end(), std::make_pair(i, j1)) == identified.end())
			identified.push_back(std::make_pair(i, j1));

		int deg_i = MG.tracker[i];

		if (collection_stage) {
			if (MG.tripleGreens.count({ deg_i,j1,j2 }) == 0 && std::find(triples_to_be_done.begin(), triples_to_be_done.end(), std::array{ deg_i,j1,j2 }) == triples_to_be_done.end()) {
				triples_to_be_done.push_back({ deg_i,j1,j2 });
			}
			return basis; //won't be used as this is collection stage, just need to return something nonempty
		}
		auto G_i_j1_j2 = MG.tripleGreens.at({ deg_i, j1, j2 });
		auto basis_product = G_i_j1_j2.getNormalBasis(0, basis);

		if (basis_product.isZero() || basis_product.size() == 1) {
			int k = std::find(distinguished.second.begin(), distinguished.second.end(), basis_product) - distinguished.second.begin();
			std::cout << distinguished.second << "\n\n";
			std::cout << "\n\n" << deg_i;
			std::cout << "\n\n" << j1;
			std::cout << "\n\n" << j2;

			return candidates[k];
		}
		auto cand_prod = id_candidates(basis_product, G_i_j1_j2.boxID, MG.NonZeroHomology.at(MG.degree[MG.index_product(deg_ij1, j2)]));

		for (auto k : distinguished.second) {
			auto it = std::find(cand_prod.begin(), cand_prod.end(), k);
			if (it != cand_prod.end())
				return candidates[it - cand_prod.begin()];
		}
		return rank_t(); //guaranteed not to happen, but will silence warnings
	}


	template<typename group_t>
	auto MultIdentify<group_t>::distinguish(int degree_ind, const std::vector<rank_t>& candidates) -> std::pair<int, std::vector<rank_t>> {

		std::vector<rank_t> products;
		products.reserve(candidates.size());
		std::vector<int> el;
		el.reserve(candidates.size());
		for (const auto& i : candidates)
			el.push_back(MG.antielement.at(std::make_pair(degree_ind, i)));
		for (int j = 0; j < MG.basicIrreducibles.size(); j++) {
			products.clear();
			auto degree_prod = MG.index_product(degree_ind, j);
			if (degree_prod == -1)
				continue;
			auto Id = MG.NonZeroHomology.at(MG.degree[degree_prod]);
			rank_t zero_matrix = rank_t::Zero(Id.group.number_of_summands());
			bool inadmissible = 0;
			for (const auto& i : el) {
				bool found = 0;
				for (auto it = MG.graph.neighbors(i).begin(); it != MG.graph.neighbors(i).end(); ++it)
					if (it.edge().id == j) { //i multiples with j to something
						products.push_back(MG.element[it.node()]); //pushback element[something]
						found = 1;
						break;
					}
				if (!found) {
					if (std::find(MG.zeroproduct[i].begin(), MG.zeroproduct[i].end(), j) == MG.zeroproduct[i].end()) {
						inadmissible = 1;
						break;
					} //else it's 0
					products.push_back(zero_matrix);
				}
			}
			if (inadmissible)
				continue;
			if (mackey::distinguish(products, Id))
				return std::make_pair(j, products);
		}
		products.clear();
		return std::make_pair(-1, products);
	}

}

