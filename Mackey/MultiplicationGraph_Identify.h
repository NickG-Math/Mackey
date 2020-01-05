#pragma once
#include "MultiplicationGraph.h"

///@file
///@brief Contains extra identification methods for the multiplication graph.

namespace Mackey {

	/// Provides extra identification methods using triple box products
	template<typename rank_t, typename diff_t>
	class MultiplicationGraphIdentify : public MultiplicationGraph<rank_t, diff_t> {
	protected:

		/// Uses triple box products to possibly identify the products of disconnected generators if identification failed before on them
		void pass_disconnected_product(bool serialize_each_step);

		/// Uses triple box products to possibly identify the products landing in degrees of disconnected generators if identification failed before on them
		void pass_disconnected_division(bool serialize_each_step);


		/// Uses triple box products to possibly identify ALL instances where identification failed before. Use with care
		void pass_all_unidentified(bool serialize_each_step);

		bool can_do_more; ///<1 if there are more triple products that can be computed for the disconnected generators
#ifdef CEREALIZE
		///Uses the Multiplication Graph constructor
		MultiplicationGraphIdentify(MultiplicationTable<rank_t, diff_t>& M) : MultiplicationGraph<rank_t, diff_t>(M) {
			triples_to_be_done.reserve(this->unidentified.size());
			identified.reserve(this->unidentified.size());
		}
#endif
		///Uses the Multiplication Graph constructor
		MultiplicationGraphIdentify(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles)
			: MultiplicationGraph<rank_t, diff_t>(level, minsphere, maxsphere, basicIrreducibles) {
			triples_to_be_done.reserve(this->unidentified.size());
			identified.reserve(this->unidentified.size());
		}

	private:
		void pass_identified(bool);
		std::vector<std::array<int, 3>> triples_to_be_done;
		///Multithreaded
		void do_triples(bool);
		std::vector<std::pair<int, int>> identified;
		void delete_identified();
		/// Uses triple box products to possibly identify given instance where identification failed before.
		bool pass_triple(std::pair<int, int>, bool);
		std::pair<int, std::vector<rank_t>> distinguish(int, const std::vector<rank_t>&);
		rank_t determine_connection_identify(const Green<rank_t, diff_t>&, int, int, int, bool);
	};


	template<typename rank_t, typename diff_t>
	void MultiplicationGraphIdentify<rank_t, diff_t>::pass_all_unidentified(bool serialize) {
		int size = 0;
		while (size != this->unidentified.size() && this->unidentified.size() > 0) {
			size = this->unidentified.size();
			for (int i = 0; i < size; i++) {
				auto thepair = this->unidentified[i];
				pass_triple(thepair, 1); //this adds extra unidentified's but no more triples so it's fine to stop at size.
			}
			pass_identified(serialize);
		}
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraphIdentify<rank_t, diff_t>::pass_disconnected_product(bool serialize) {
		can_do_more = 0;
		for (const auto& i : this->disconnected) {
			for (int j = 0; j < this->number_of_irreducibles; j++) {
				auto ij = std::make_pair(i, j);
				if (find(this->unidentified, ij) != -1 && pass_triple(ij, 1)) { //we may be able to identify using one j so break
					if (j != this->number_of_irreducibles - 1)
						can_do_more = 1; //one per j might suffice, if it doesn't can_do_more let's us know that we can go again
					break;
				}
			}
		}
		pass_identified(serialize);
	}


	template<typename rank_t, typename diff_t>
	void MultiplicationGraphIdentify<rank_t, diff_t>::pass_disconnected_division(bool serialize) {
		can_do_more = 1;
		for (const auto& pair : this->unidentified) {
			auto deg = this->index_product(this->tracker[pair.first], pair.second);
			for (const auto& i : this->disconnected) {
				if (this->tracker[i] == deg && pass_triple(pair, 1))
					break;
			}
			can_do_more = 0;
		}
		pass_identified(serialize);
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraphIdentify<rank_t, diff_t>::pass_identified(bool serialize) {
		if (triples_to_be_done.size() != 0)
			do_triples(serialize);
		for (const auto& pair : identified)
			pass_triple(pair, 0);
		delete_identified();
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraphIdentify<rank_t, diff_t>::delete_identified() {
		std::vector<std::pair<int, int>> true_unidentified;
		true_unidentified.reserve(this->unidentified.size() - identified.size());
		for (const auto& i : this->unidentified) {
			if (find(identified, i) == -1)
				true_unidentified.push_back(i);
		}
		this->unidentified = true_unidentified;
		identified.clear();
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraphIdentify<rank_t, diff_t>::do_triples(bool serialize) {
#pragma omp parallel for num_threads(12) schedule(dynamic)
		for (int i = 0; i < triples_to_be_done.size(); i++) {
			this->triple_product(triples_to_be_done[i][0], triples_to_be_done[i][1], triples_to_be_done[i][2]);
		}
		triples_to_be_done.clear();
#ifdef CEREALIZE
		if (serialize)
			saver(*this, "Table_with_triples.bin", "binary");
#endif
	}

	template<typename rank_t, typename diff_t>
	bool MultiplicationGraphIdentify<rank_t, diff_t>::pass_triple(std::pair<int, int> ij, bool collection_stage) {
		auto i = ij.first;
		auto j = ij.second;
		auto deg_i = this->tracker[i];
		auto G = this->Greens[deg_i][j];
		int deg_ij = this->index_product(deg_i, j); //can't be -1
		auto size_now = identified.size();
		auto v = determine_connection_identify(G, i, j, deg_ij, collection_stage);

		if (collection_stage && v.size() != 0)
			return 1;
		else if (v.size() != 0) {
			this->connect(v, i, j, deg_ij);
			return 1;
		}
		return 0;
	}

	template<typename rank_t, typename diff_t>
	rank_t MultiplicationGraphIdentify<rank_t, diff_t>::determine_connection_identify(const Green<rank_t, diff_t>& G_ij1, int i, int j1, int deg_ij1, bool collection_stage) {
		rank_t basis = G_ij1.getNormalBasis(0, this->element[i]);

		std::vector<rank_t> candidates = id_candidates(basis, G_ij1.boxID, this->NonZeroHomology[deg_ij1]);

		for (const auto& cand : candidates) { //add all candidates and their edges
			int t= this->find_and_add_element(deg_ij1, cand);
			if (t != this->element.size() - 1) //nothing added
				continue;
			this->add_edges(t);
		}

		std::pair<int, std::vector<rank_t>> distinguished = distinguish(deg_ij1, candidates); //see if candidates can be distinguished from the edges out of them
		int j2 = distinguished.first;
		if (j2 == -1)
			return rank_t();

		if (find(identified, std::make_pair(i, j1))== -1)
			identified.push_back(std::make_pair(i, j1));

		int deg_i = this->tracker[i];

		if (collection_stage) {
			if (this->tripleGreens.count({ deg_i,j1,j2 }) == 0 && find<std::vector<std::array<int, 3>>, std::array<int, 3>>(triples_to_be_done, { deg_i,j1,j2 }) == -1)
				triples_to_be_done.push_back({ deg_i,j1,j2 });
			return basis; //won't be used as this is collection stage, just need to return something nonempty
		}

		auto G_i_j1_j2 = this->tripleGreens.at({ deg_i, j1, j2 });
		auto basis_product = G_i_j1_j2.getNormalBasis(0, basis);

		if (basis_product.isZero() || basis_product.size() == 1) {
			int k = find(distinguished.second, basis_product);
			return candidates[k];
		}
		auto cand_prod = id_candidates(basis_product, G_i_j1_j2.boxID, this->NonZeroHomology[this->index_product(deg_ij1, j2)]);

		for (int k = 0; k < distinguished.second.size(); k++) {
			auto s = find(cand_prod, distinguished.second[k]);
			if (s != -1)
				return candidates[s];
		}
		return rank_t(); //guaranteed not to happen, but will silence warnings
	}


	template<typename rank_t, typename diff_t>
	std::pair<int, std::vector<rank_t>> MultiplicationGraphIdentify<rank_t, diff_t>::distinguish(int degree_ind, const std::vector<rank_t>& candidates) {

		std::vector<rank_t> products;
		products.reserve(candidates.size());
		std::vector<int> el;
		el.reserve(candidates.size());
		for (const auto& i : candidates) {
			el.push_back(this->antielement.at(std::make_pair(degree_ind, i)));
		}
		for (int j = 0; j < this->number_of_irreducibles; j++) {
			products.clear();
			auto degree_prod = this->index_product(degree_ind, j);
			if (degree_prod == -1)
				continue;
			auto Id = this->NonZeroHomology[degree_prod];
			rank_t zero_matrix = rank_t::Zero(Id.group.size());
			bool inadmissible = 0;
			for (const auto& i : el) {
				auto k = find(this->edgeid[i], j);
				if (k == -1) {
					if (find(this->zeroproduct[i], j) == -1) {
						inadmissible = 1;
						break;
					} //else it's 0
					products.push_back(zero_matrix);
				}
				else
					products.push_back(this->element[this->edges[i][k]]);
			}
			if (inadmissible)
				continue;
			if (Mackey::distinguish(products, Id))
				return std::make_pair(j, products);
		}
		products.clear();
		return std::make_pair(-1, products);
	}

}

