#pragma once
#include "MultiplicationGraph.h"

///@file
///@brief Contains extra identification methods for the multiplication graph.

namespace Mackey {

	/// Provides extra identification methods using triple box products
	template<typename rank_t, typename diff_t>
	class MultiplicationGraphConnectivity : public MultiplicationGraph<rank_t, diff_t> {
	public:
#ifdef CEREALIZE
		///Uses the Multiplication Graph constructor
		MultiplicationGraphConnectivity(MultiplicationTable<rank_t, diff_t>& M) : MultiplicationGraph<rank_t, diff_t>(M) {}
#endif
		///Uses the Multiplication Graph constructor
		MultiplicationGraphConnectivity(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles)
			: MultiplicationGraph<rank_t, diff_t>(level, minsphere, maxsphere, basicIrreducibles) {}

	///The generators disconnected from the sources (element indices)
	std::vector<int> trully_disconnected;

	///The degrees of the trully_disconnected
	std::vector<std::vector<int>> trully_disconnected_degrees;

	///Computes the nodes of the Multiplication Graph disconnected from the given sources
	void compute_with_sources(const std::vector<std::vector<int>>&);

	private:
		Graph<void> Wrong; ///has same connectivity as multiplication graph but potentially "wrong" edges 
		std::vector<int> sources;
		void set_sources(const std::vector<std::vector<int>>&);
		void add_candidates(const Green<rank_t, diff_t>&, int, int, int);
		void force_connection(const Green<rank_t, diff_t>&, int, int, int);
	};


	template<typename rank_t, typename diff_t>
	void MultiplicationGraphConnectivity<rank_t, diff_t>::compute_with_sources(const std::vector<std::vector<int>>& given_sources) {
		set_sources(given_sources);
		for (const auto& i : sources)
			this->computeWithSource(i);
		this->compute_disconnected();
		for (const auto& i : this->disconnected) {
			for (int j = 0; j < this->number_of_irreducibles; j++) {
				auto ij = std::make_pair(i, j);
				if (find(this->unidentified, ij) != -1) {
					auto deg_i = this->tracker[i];
					auto G = this->Greens[deg_i][j];
					auto deg_ij = this->index_product(deg_i, j);
					add_candidates(G, i, j, deg_ij);
				}
			}
		}
		Wrong.number_of_nodes = this->number_of_nodes;
		Wrong.edges.resize(this->number_of_nodes);
		for (int i = 0; i < this->number_of_nodes; i++)
			Wrong.edges[i] = this->edges[i];

		for (const auto& i : this->disconnected) {
			for (int j = 0; j < this->number_of_irreducibles; j++) {
				auto ij = std::make_pair(i, j);
				if (find(this->unidentified, ij) != -1) {
					auto deg_i = this->tracker[i];
					auto G = this->Greens[deg_i][j];
					auto deg_ij = this->index_product(deg_i, j);
					force_connection(G, i, j, deg_ij);
				}
			}
		}

		for (const auto& i : sources)
			Wrong.computeWithSource(i);

		Wrong.compute_disconnected(this->number_of_generators);
		trully_disconnected = Wrong.disconnected;

		trully_disconnected_degrees.reserve(trully_disconnected.size());
		for (const auto& i : trully_disconnected)
			trully_disconnected_degrees.push_back(this->getdegree(i));
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraphConnectivity<rank_t, diff_t>::set_sources(const std::vector<std::vector<int>>& given_sources) {
		sources.reserve(given_sources.size());
		for (const auto& i : given_sources) {
			auto deg = this->antidegree[i];
			auto basis = basisElement<rank_t>(1, 0);
			sources.push_back(this->antielement[std::make_pair(deg, basis)]);
		}
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraphConnectivity<rank_t, diff_t>::add_candidates(const Green<rank_t, diff_t>& G_ij, int i, int j, int deg_ij) {
		rank_t basis = G_ij.getNormalBasis(0, this->element[i]);
		std::vector<rank_t> candidates = id_candidates(basis, G_ij.boxID, this->NonZeroHomology[deg_ij]);
		for (const auto& cand : candidates) { //add all candidates and their edges
			int t = this->find_and_add_element(deg_ij, cand);
			if (t != this->element.size() - 1) //nothing added
				continue;
			this->add_edges(t);
		}
	}



	template<typename rank_t, typename diff_t>
	void MultiplicationGraphConnectivity<rank_t, diff_t>::force_connection(const Green<rank_t, diff_t>& G_ij, int i, int j, int deg_ij) {
		rank_t basis = G_ij.getNormalBasis(0, this->element[i]);
		std::vector<rank_t> candidates = id_candidates(basis, G_ij.boxID, this->NonZeroHomology[deg_ij]);
		std::vector<rank_t> notconnected;
		notconnected.reserve(candidates.size());
		for (const auto& cand : candidates) {
			int t = this->antielement.at(std::make_pair(deg_ij, cand));
			//check if candidate is connected
			if (find(this->disconnected, t) != -1)
				notconnected.push_back(cand);
			//check if we have an injection
			if (!this->injection(i, t))
				return;
		}
		if (!notconnected.empty()) {
			//get all connected elements
			std::vector<rank_t> connected;
			connected.reserve(candidates.size());
			for (int t = 0; t < this->number_of_nodes; t++) {
				if (this->tracker[t] == deg_ij && find(this->disconnected, t) == -1)
					connected.push_back(this->element[t]);
			}
			//check if the non connected ones are in the span of the connected ones, hence connected
			for (const auto& cand : notconnected) {
				if (!inSpan(cand, connected, this->NonZeroHomology[deg_ij].group))
					return; //if not then don't bother
			}
			//connect the first connected and first candidate
			Wrong.edges[this->antielement.at(std::make_pair(deg_ij, connected[0]))].push_back(this->antielement.at(std::make_pair(deg_ij, candidates[0])));
		}
		//connect first candidates to i
		Wrong.edges[this->antielement.at(std::make_pair(deg_ij, candidates[0]))].push_back(i);
	}
}
