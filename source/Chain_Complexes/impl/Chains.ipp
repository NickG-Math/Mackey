#pragma once
#include "../Chains.hpp"

///@file
///@brief Contains Chains and Junction classes for using chain complexes.

namespace mackey {


	template<typename rank_t, typename diff_t>
	Arrow<rank_t, diff_t>::Arrow(const rank_t& domain, const rank_t& range, const diff_t& diff) : domain(domain), range(range), diff(diff) {}


	//The maximum index of the Chain complex
	template<typename rank_t, typename diff_t>
	int Chains<rank_t, diff_t>::maxindex() const {
		return (int)rank.size() - 1;
	}

	template<typename rank_t, typename diff_t>
	Chains<rank_t, diff_t>::Chains(const std::vector<rank_t>& rank, const std::vector<diff_t>& diff) : rank(rank), diff(diff) {}	///<Constructor given ranks and diffs

	template<typename rank_t, typename diff_t>
	void Chains<rank_t, diff_t>::reserve(size_t i) {
		rank.reserve(i);
		diff.reserve(i);
	}

	template<typename rank_t, typename diff_t>
	void Chains<rank_t, diff_t>::push_back(const Arrow<rank_t, diff_t>& arrow) {
		rank.push_back(arrow.domain);
		diff.push_back(arrow.diff);
	}

	template<typename rank_t, typename diff_t>
	void Chains<rank_t, diff_t>::push_back(const rank_t& domain, const diff_t& differential) {
		rank.push_back(domain);
		diff.push_back(differential);
	}

	//The dual cochain complex reindexed as a Chain complex. Optionally stop at index k (if k=-1 then nonstop i.e. k=maxindex)
	template<typename rank_t, typename diff_t>
	Chains<rank_t, diff_t> Chains<rank_t, diff_t>::dualize(int k) const {
		if (k == -1)
			k = maxindex();
		Chains dual;
		dual.reserve(maxindex() + 1);
		dual.push_back(rank[maxindex()], diff_t());
		for (int i = 1; i <= k; i++)
			dual.push_back(rank[maxindex() - i], diff[maxindex() - i + 1].transpose());
		return dual;
	}

	template<typename rank_t, typename diff_t>
	bool Chains<rank_t, diff_t>::operator==(const Chains& C) const {
		return (rank == C.rank && diff = C.diff);
	}

	template<typename rank_t, typename diff_t>
	std::ostream& operator<< (std::ostream& os, const Chains<rank_t, diff_t>& C) {
		for (size_t i = 1; i < C.diff.size(); i++) {
			os << "diff[" << i << "]: [";
			os << (C.rank[i].template cast<int64_t>()) << "] -> [";
			os << C.rank[i - 1] .template cast<int64_t>() << "] is:\n";
			os << C.diff[i].template cast<int64_t>() << "\n";
		}
		return os;
	}


	//Extract the Junction C[i+1]->C[i]->C[i-1] from the Chains C.
	template<typename rank_t, typename diff_t>
	Junction<rank_t, diff_t>::Junction(const Chains<rank_t, diff_t>& C, int i) {
		rank = C.rank[i];
		if (0 < i)
			setOut(C.rank[i - 1], C.diff[i]);
		if (i < C.maxindex())
			setIn(C.rank[i + 1], C.diff[i + 1]);
	}
	template<typename rank_t, typename diff_t>
	Junction<rank_t, diff_t>::Junction(const Arrow<rank_t, diff_t>& In, const Arrow<rank_t, diff_t>& Out) {
		setIn(In);
		setOut(Out, 0);
	}

	//Construct Junction by directly setting the elements
	template<typename rank_t, typename diff_t>
	Junction<rank_t, diff_t>::Junction(const rank_t& rank, const rank_t& rankOut, const rank_t& rankIn, const diff_t& diffOut, const diff_t& diffIn)
		: rank(rank), rankOut(rankOut), rankIn(rankIn), diffOut(diffOut), diffIn(diffIn) {}

	template<typename rank_t, typename diff_t>
	void Junction<rank_t, diff_t>::setIn(const Arrow<rank_t, diff_t>& arrow, bool docenter) {
		rankIn = arrow.domain;
		diffIn = arrow.diff;
		if (docenter)
			rank = arrow.range;
	}

	template<typename rank_t, typename diff_t>
	void Junction<rank_t, diff_t>::setIn(const rank_t& domain, const diff_t& diff) {
		rankIn = domain;
		diffIn = diff;
	}

	template<typename rank_t, typename diff_t>
	void Junction<rank_t, diff_t>::setOut(const Arrow<rank_t, diff_t>& arrow, bool docenter) {
		rankOut = arrow.range;
		diffOut = arrow.diff;
		if (docenter)
			rank = arrow.domain;
	}

	template<typename rank_t, typename diff_t>
	void Junction<rank_t, diff_t>::setOut(const rank_t& range, const diff_t& diff) {
		rankOut = range;
		diffOut = diff;
	}

}
