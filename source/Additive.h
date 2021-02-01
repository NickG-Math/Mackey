#pragma once
#include "Compute.h"

///@file
///@brief Contains the class for presenting the additive structure of the RO(G) homology. 

namespace Mackey {

	///The additive structure of the RO(G) homology of a point
	template<typename rank_t, typename diff_t>
	class AdditiveStructure {
	public:
		std::vector<int> minsphere;///<The lower bound on the range of our spheres
		std::vector<int> maxsphere;///<The upper bound on the range of our spheres

		///Construct given the range of our spheres
		AdditiveStructure(const std::vector<int>& minsphere, const std::vector<int>& maxsphere) :minsphere(minsphere), maxsphere(maxsphere) {
			compute();
		}

#ifdef CEREALIZE
		///Default constructor for serealization
		AdditiveStructure() {}
		template<typename Archive, typename s_rank_t, typename s_diff_t>
		friend void serialize(Archive&, AdditiveStructure<s_rank_t, s_diff_t>&);
#endif

		///Retrieve the Mackey functors in the homology of the given sphere
		template<typename sphere_t>
		std::vector<MackeyFunctor<rank_t>> getMackey(const sphere_t&) const;

		///Retrieve the Mackey functor for the given degree
		template<typename sphere_t>
		MackeyFunctor<rank_t> getMackey(int, const sphere_t&) const;


		///Print the Mackey functors for each degree
		template<typename T>
		void print_answer(T& stream);

		///Print the unique Mackey functors that appear
		template<typename T>
		void print_unique(T& stream);

		///Print the unnamed Mackey functors that appear
		template<typename T>
		void print_unknown(T& stream);

		///Identify the Mackey functors
		void identify();
	private:
		std::vector<int> true_unknowns;
		std::vector<MackeyFunctor<rank_t>> uniqueMackeys;
		std::vector<MackeyFunctor<rank_t>> unknownMackeys;
		std::vector<std::vector<MackeyFunctor<rank_t>>> allMackeys;
		std::vector<char> found_unknown;
		template<typename sphere_t>
		int hash(const sphere_t& sphere) const { return hashvector(sphere, minsphere, maxsphere); }
		void compute();
		void identify_unique(MackeyFunctor<rank_t>&, bool);
		void identify_unknown(int);

		void check_unique_and_assign(const MackeyFunctor<rank_t>&);
		void check_unknown_and_assign(const std::vector<MackeyFunctor<rank_t>>&);

	};

	template<typename rank_t, typename diff_t>
	void AdditiveStructure<rank_t, diff_t>::compute() {

		auto spheres = DegreeConstruction(minsphere, maxsphere);
		allMackeys.resize(hash(maxsphere) + 1);
#pragma omp parallel for num_threads(12)
		for (int i = 0; i < spheres.size(); i++)
		{
			auto M = ROHomology<rank_t, diff_t>(spheres[i]);
			int hasher = hashvector(spheres[i], minsphere, maxsphere);
			allMackeys[hasher].reserve(M.size());
			for (const auto& Mack : M)
			{
				allMackeys[hasher].push_back(Mack);
			}
		}
	}

	template<typename rank_t, typename diff_t>
	template<typename sphere_t>
	std::vector<MackeyFunctor<rank_t>> AdditiveStructure<rank_t, diff_t>::getMackey(const sphere_t& sphere) const {
		return allMackeys[hash(sphere)];
	}

	template<typename rank_t, typename diff_t>
	template<typename sphere_t>
	MackeyFunctor<rank_t> AdditiveStructure<rank_t, diff_t>::getMackey(int degree, const sphere_t& sphere) const {
		return getMackey(sphere)[Reindex(degree, sphere)];
	}


	template<typename rank_t, typename diff_t>
	template<typename T>
	void AdditiveStructure<rank_t, diff_t>::print_answer(T& stream) {
		auto spheres = DegreeConstruction(minsphere, maxsphere);
		for (const auto& i : spheres) {
			int k = 0;
			for (const auto& Mack : allMackeys[hash(i)]) {
				stream << "The k=" << invReindex(k, i) << " homology of the " << i << " sphere is " << Mack << "\n";
				k++;
			}
		}
	}

	template<typename rank_t, typename diff_t>
	template<typename T>
	void AdditiveStructure<rank_t, diff_t>::print_unique(T& stream) {
		for (const auto& i : uniqueMackeys)
			stream << i << "\n";
	}

	template<typename rank_t, typename diff_t>
	template<typename T>
	void AdditiveStructure<rank_t, diff_t>::print_unknown(T& stream) {
		for (const auto& i : true_unknowns)
			stream << unknownMackeys[i] << ": \n" << unknownMackeys[i].print() << "\n";
	}


	template<typename rank_t, typename diff_t>
	void AdditiveStructure<rank_t, diff_t>::check_unique_and_assign(const MackeyFunctor<rank_t>& M) {
		for (const auto& i : uniqueMackeys) {
			if (i.name == M.name)
				return;
		}
		uniqueMackeys.push_back(M);
	}

	template<typename rank_t, typename diff_t>
	void AdditiveStructure<rank_t, diff_t>::check_unknown_and_assign(const std::vector<MackeyFunctor<rank_t>>& isos) {
		for (const auto& i : unknownMackeys) {
			if (i.isomorphic(isos))
				return;
		}
		unknownMackeys.push_back(isos[0]);
	}


	template<typename rank_t, typename diff_t>
	void AdditiveStructure<rank_t, diff_t>::identify_unique(MackeyFunctor<rank_t>& M, bool check_unknown) {
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
		for (auto & i:uniqueMackeys) {
			for (auto& j : uniqueMackeys) {
				auto k=i+j;
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


	template<typename rank_t, typename diff_t>
	void AdditiveStructure<rank_t, diff_t>::identify_unknown(int index) {
		if (!unknownMackeys[index].name.empty())
			return;
		identify_unique(unknownMackeys[index], 0);
		if (!unknownMackeys[index].name.empty())
			return;
		auto isos = unknownMackeys[index].isomorphism_class();
		for (const auto& i : uniqueMackeys) {
			for (int j = 0; j < unknownMackeys.size(); j++) {
				if (j != index && (i + unknownMackeys[j]).isomorphic(isos)) {
					unknownMackeys[index].name = i.name + "+ unknown " + std::to_string(j);
					return;
				}
			}
		}
		for (int i = 0; i < unknownMackeys.size(); i++) {
			if (i != index) {
				for (int j = 0; j < unknownMackeys.size(); j++) {
					if (j != index && (unknownMackeys[i] + unknownMackeys[j]).isomorphic(isos)) {
						unknownMackeys[index].name = "unknown " + std::to_string(i) + "+ unknown " + std::to_string(j);
						return;
					}
				}
			}
		}
	}

	template<typename rank_t, typename diff_t>
	void AdditiveStructure<rank_t, diff_t>::identify() {
		uniqueMackeys.reserve(allMackeys.size());
		unknownMackeys.reserve(allMackeys.size());
		for (auto& i : allMackeys) {
			for (auto& j : i)
				identify_unique(j, 1);
		}
		int size;
		do {
			size = uniqueMackeys.size();
			for (int i = 0; i < unknownMackeys.size(); i++) {
				identify_unknown(i);
			}
		} while (size != uniqueMackeys.size());
		true_unknowns.reserve(unknownMackeys.size());
		for (int i = 0; i < unknownMackeys.size(); i++) {
			if (unknownMackeys[i].name.empty()) {
				unknownMackeys[i].name = "unknown " + std::to_string(i);
				true_unknowns.push_back(i);
			}
		}
		for (auto& i : allMackeys) {
			for (auto& j : i) {
				if (j.name.empty()) {
					for (const auto& k : unknownMackeys) {
						if (j.isomorphic(k))
							j.name = k.name;
					}
				}
			}
		}
	}

}


