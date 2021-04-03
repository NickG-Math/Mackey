#pragma once
#include "../Mult_Connectivity.hpp"

///@file
///@brief Contains methods for the finding the connectivity of the multiplication graph should identification fail.

namespace mackey {

	template<typename group_t>
	template<typename...Args>
	MultConnectivity<group_t>::MultConnectivity(Args&&... args) : MultGraph<group_t>(std::forward<Args>(args)...), shortest_paths(this->graph) {}

	template<typename group_t>
	void MultConnectivity<group_t>::set_sources(const std::vector<std::vector<int>>& given_sources) {
		sources.reserve(given_sources.size());
		for (const auto& i : given_sources) {
			auto deg = this->antidegree[i];
			rank_t basis(1);
			basis[0] = 1;
			sources.push_back(this->antielement.at(std::make_pair(deg, basis)));
		}
	}

	template<typename group_t>
	void MultConnectivity<group_t>::compute_with_sources(const std::vector<std::vector<int>>& given_sources) {
		set_sources(given_sources);
		for (const auto& i : sources)
			shortest_paths.compute_with_root(i);
		for (int i = 0; i < this->graph.number_of_nodes(); i++)
			if (shortest_paths.disconnected(i))
				before_disconnected.insert(i);
		doinstages(1); //first collection stage
		doinstages(0); //the force connection stage
		for (const auto& i : sources)
			shortest_paths.compute_with_root(i);
		for (int i = 0; i < this->number_of_generators; i++) {
			if (shortest_paths.disconnected(i)) {
				disconnected_indices.push_back(i);
				disconnected_degrees.push_back(this->getdegree(i));
			}
		}
	}

	template<typename group_t>
	void MultConnectivity<group_t>::doinstages(bool stage) {
		for (auto i : before_disconnected)
			for (int j = 0; j < this->basicIrreducibles.size(); j++) 
				if (this->unidentified.find(std::pair(i,j))!=this->unidentified.end()) 
					wrong_connection(i, j, stage);
	}

	template<typename group_t>
	void MultConnectivity<group_t>::wrong_connection(int i, int j, bool collection_stage) {
		auto prod=this->product_element_irreducible(i, j);
		auto candidates = this->product_candidates(prod);
		if (collection_stage) {
			for (const auto& cand : candidates) { //add all candidates and their edges
				int t = this->find_and_add_element(prod.deg, cand);
				if (t != this->element.size() - 1) //nothing added
					continue;
				this->add_edges(t);
			}
			return;
		}
		std::vector<rank_t> notconnected;
		notconnected.reserve(candidates.size());
		for (const auto& cand : candidates) {
			int t = this->antielement.at(std::make_pair(prod.deg, cand));
			//check if candidate is connected
			if (before_disconnected.find(t) != before_disconnected.end())
				notconnected.push_back(cand);
			//check if we have an injection
			if (!this->injection(i, t))
				return;
		}
		if (!notconnected.empty()) {
			//get all connected elements
			std::vector<rank_t> connected;
			connected.reserve(candidates.size());
			for (int t = 0; t < this->graph.number_of_nodes(); t++) {
				if (this->tracker[t] == prod.deg && before_disconnected.find(t) == before_disconnected.end())
					connected.push_back(this->element[t]);
			}
			//check if the non connected ones are in the span of the connected ones, hence connected
			for (const auto& cand : notconnected) {
				if (this->NonZeroHomology.at(this->degree[prod.deg]).group.span(cand,connected).size()==0)
					return; //if it's not in the span then don't bother
			}
			//connect the first connected to first candidate
			auto cond_ind = this->antielement.at(std::make_pair(prod.deg, connected[0]));
			auto cand_ind = this->antielement.at(std::make_pair(prod.deg, candidates[0]));
			this->graph.draw_edge(cond_ind, cand_ind, 0, j);
		}
		//connect first candidate to i
		auto cand_ind = this->antielement.at(std::make_pair(prod.deg, candidates[0]));
		this->graph.draw_edge(cand_ind, i, 0, j);
	}
}
