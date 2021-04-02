#pragma once
#include "../Additive.hpp"

///@file
///@brief Contains the class for presenting the additive structure of the RO(G) homology. 

namespace mackey {

	///Construct given the range of our spheres
	template<typename group_t>
	AdditiveStructure<group_t>::AdditiveStructure(const std::vector<int>& minsphere, const std::vector<int>& maxsphere) :minsphere(minsphere), maxsphere(maxsphere) {
		compute();
		identify();
	}

	///Retrieve the Mackey functors in the homology of the given sphere
	template<typename group_t>
	template<typename sphere_t>
	std::vector<MackeyFunctor<typename group_t::rank_t>> AdditiveStructure<group_t>::getMackey(const sphere_t& sphere) const {
		return allMackeys[perfect_hash(sphere, minsphere, maxsphere)];
	}

	///Retrieve the Mackey functor for the given degree
	template<typename group_t>
	template<typename sphere_t>
	MackeyFunctor<typename group_t::rank_t> AdditiveStructure<group_t>::getMackey(int degree, const sphere_t& sphere) const {
		return getMackey(sphere)[Reindex<group_t>(degree, sphere)];
	}

	template<typename group_t>
	void AdditiveStructure<group_t>::compute() {
		auto spheres = InterpolatingVectorGenerator<std::vector<int>>(minsphere, maxsphere)();
		allMackeys.resize(perfect_hash(maxsphere, minsphere, maxsphere) + 1);
		MACKEY_RUN_LOOP_PARALLEL
			for (int i = 0; i < spheres.size(); i++)
				allMackeys[perfect_hash(spheres[i], minsphere, maxsphere)] = ROHomology<group_t>(spheres[i]);
	}


	template<typename group_t>
	std::ostream& operator<<(std::ostream& os, const AdditiveStructure<group_t>& A) {
		InterpolatingVectorGenerator<std::vector<int>> spheres(A.minsphere, A.maxsphere);
		for (const auto& i : spheres) {
			int k = 0;
			for (const auto& Mack : A.allMackeys[perfect_hash(i, A.minsphere, A.maxsphere)]) {
				std::vector<int> deg;
				deg.push_back(k);
				deg.insert(deg.end(),i.begin(), i.end());
				os << "The k=" << std::setw(4)<< invReindex<group_t>(deg).front() << "  homology of the  " << std::setw(2) << i << " sphere is " << Mack << "\n";
				k++;
			}
		}
		return os;
	}

	template<typename group_t>
	const std::vector<MackeyFunctor<typename group_t::rank_t>>& AdditiveStructure<group_t>::identified() {
		return uniqueMackeys;
	}


	template<typename group_t>
	const std::vector<MackeyFunctor<typename group_t::rank_t>>& AdditiveStructure<group_t>::unknown() {
		return unknownMackeys;
	}

	template<typename group_t>
	void AdditiveStructure<group_t>::check_unique_and_assign(const MackeyFunctor<rank_t>& M) {
		for (const auto& i : uniqueMackeys) 
			if (i.name == M.name)
				return;
		uniqueMackeys.push_back(M);
	}

	template<typename group_t>
	void AdditiveStructure<group_t>::check_unknown_and_assign(const std::vector<MackeyFunctor<rank_t>>& isos) {
		for (const auto& i : unknownMackeys) 
			if (i.isomorphic(isos))
				return;
		unknownMackeys.push_back(isos[0]);
	}


	template<typename group_t>
	void AdditiveStructure<group_t>::identify_unique(MackeyFunctor<rank_t>& M, bool check_unknown) {
		M.notation();
		if (!M.name.empty()) {
			check_unique_and_assign(M);
			return;
		}
		auto isos = M.isomorphism_class();
		for (const auto& i : uniqueMackeys) {
			if (i.isomorphic(isos)) {
				M.name = i.name;
				return;
			}
		}
		for (auto& i : uniqueMackeys) {
			for (auto& j : uniqueMackeys) {
				auto k = i + j;
				if (k.isomorphic(isos)) {
					M.name = k.name;
					uniqueMackeys.push_back(M);
					return;
				}
			}
		}
		if (check_unknown)
			check_unknown_and_assign(isos);
	}


	template<typename group_t>
	void AdditiveStructure<group_t>::identify_unknown(int index) {
		if (!unknownMackeys[index].name.empty())
			return;
		identify_unique(unknownMackeys[index], 0);
		if (!unknownMackeys[index].name.empty())
			return;
		auto isos = unknownMackeys[index].isomorphism_class();
		for (const auto& i : uniqueMackeys) {
			for (size_t j = 0; j < unknownMackeys.size(); j++) {
				if (j != index && (i + unknownMackeys[j]).isomorphic(isos)) {
					unknownMackeys[index].name = i.name + "+ unknown " + std::to_string(j);
					return;
				}
			}
		}
		for (size_t i = 0; i < unknownMackeys.size(); i++)
			if (i != index)
				for (size_t j = 0; j < unknownMackeys.size(); j++)
					if (j != index && (unknownMackeys[i] + unknownMackeys[j]).isomorphic(isos)) {
						unknownMackeys[index].name = "unknown " + std::to_string(i) + "+ unknown " + std::to_string(j);
						return;
					}
	}

	template<typename group_t>
	void AdditiveStructure<group_t>::identify() {
		uniqueMackeys.reserve(allMackeys.size());
		unknownMackeys.reserve(allMackeys.size());
		for (auto& i : allMackeys) {
			for (auto& j : i)
				identify_unique(j, 1);
		}
		int size;
		do {
			size = uniqueMackeys.size();
			for (size_t i = 0; i < unknownMackeys.size(); i++)
				identify_unknown(i);
		} while (size != uniqueMackeys.size());
		true_unknowns.reserve(unknownMackeys.size());
		for (size_t i = 0; i < unknownMackeys.size(); i++) {
			if (unknownMackeys[i].name.empty()) {
				unknownMackeys[i].name = "unknown " + std::to_string(i);
				true_unknowns.push_back(i);
			}
		}
		for (auto& i : allMackeys)
			for (auto& j : i)
				if (j.name.empty())
					for (const auto& k : unknownMackeys)
						if (j.isomorphic(k))
							j.name = k.name;
	}

}


