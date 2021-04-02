#pragma once
#include "../Mult_Table.hpp"

///@file
///@brief Contains the class and methods to form the multiplication table.

namespace mackey
{
	template<typename group_t>
	MultTable<group_t>::MultTable(int _level, const std::vector<int>& _minsphere, const std::vector<int>& _maxsphere, const std::vector<std::vector<int>>& _basicIrreducibles)
		: perfect_hasher(_minsphere, _maxsphere)  {
		level = _level;
		minsphere = _minsphere;
		maxsphere = _maxsphere;
		basicIrreducibles = _basicIrreducibles;
		getIrreducibleChains();
		make();
	}
	
	template<typename group_t>
	MultTable<group_t>::MultTable(const MultTableData<group_t>& MTD) : MultTableData<group_t>(MTD), perfect_hasher(minsphere, maxsphere) {
		getIrreducibleChains();
	}

	template<typename group_t>
	MultTable<group_t>::MultTable(MultTableData<group_t>&& MTD) : MultTableData<group_t>(std::move(MTD)), perfect_hasher(minsphere, maxsphere) {
		getIrreducibleChains();
	}


	template<typename group_t>
	int MultTable<group_t>::getdegreeindex(const std::vector<int>& degree) const {
		auto iterator = this->antidegree.find(degree);
		if (iterator == this->antidegree.end())
			return -1;
		else
			return iterator->second;
	}

	template<typename group_t>
	bool MultTable<group_t>::degreewithinrange(const std::vector<int> &degree) const {
		for (size_t i = 1; i < degree.size(); i++)
			if (degree[i]<minsphere[i - 1] || degree[i]>maxsphere[i-1])
				return 0;
		return 1;
	}

	///Set the basicChains given the basicIrreducibles
	template<typename group_t>
	void MultTable<group_t>::getIrreducibleChains() {
		basicChains.reserve(basicIrreducibles.size());
		for (const auto& i : basicIrreducibles)
			basicChains.push_back(ROChains<group_t>(getsphere(i)));
	}


	/// The Multiplication Table
	template<typename group_t>
	std::vector<int> MultTable<group_t>::getsphere(const std::vector<int>& deg) const {
		return tail(deg.data(), deg.size(), 1);
	}

	///Set the basicChains given the basicIrreducibles
	template<typename group_t>
	void MultTable<group_t>::generators() {

		InterpolatingVectorGenerator<std::vector<int>> spheres(minsphere, maxsphere);
		auto sphere_it = spheres.begin();
		MACKEY_RUN_LOOP_PARALLEL
		for (int k = 0; k < spheres.size(); k++) {
			std::vector<int> sphere;
			MACKEY_RUN_BLOCK_SERIAL{
				sphere = *sphere_it;
				++sphere_it;
			}
			auto C = ROChains<group_t>(sphere);
			std::vector<int> deg;
			deg.reserve(1 + sphere.size());
			for (int k = 0; k <= C.maxindex();k++) {
				Junction<rank_t, diff_t> J(C, k);
				internal::IDGeneratorCompute<group_t>ID(level, J);
				if (ID.ID.group.number_of_summands() != 0) {
					deg.clear();
					deg.push_back(k);
					for (const auto& j : sphere)
						deg.push_back(j);
					deg = invReindex<group_t>(deg);
					MACKEY_RUN_BLOCK_SERIAL{
						NonZeroHomology[deg] = ID.ID;
					}
				}
			}
		}
		//reinterpret map into vector
		degree.reserve(NonZeroHomology.size());
		int degindex = 0;
		for (auto it = NonZeroHomology.begin(); it != NonZeroHomology.end(); it++) {
			auto deg = it->first;
			degree.push_back(deg);
			antidegree[deg] = degindex;
			degindex++;
		}

	}

	///Compute the products of all elements and irreducibles
	template<typename group_t>
	void MultTable<group_t>::make() {
		generators();
		multiply_all_indices();

		std::vector<std::pair<int, int>> within_range;
		within_range.reserve(index_product.size());
		for (size_t j = 0; j < basicIrreducibles.size(); j++) {
			for (size_t i = 0; i < degree.size(); i++)
				if (index_product(i, j) != -1)
					within_range.push_back(std::make_pair(i, j));
		}
		compute(within_range);
	}



	template<typename group_t>
	void MultTable<group_t>::compute(const std::vector<std::pair<int, int>>& pairs) {
		for (const auto& i : pairs) {
			auto sphere = getsphere(degree[i.first]);
			auto j = perfect_hasher(sphere);
			if (sphereChains.count(j) == 0)
				sphereChains[j] = ROChains<group_t>(sphere);
		}
		MACKEY_RUN_LOOP_PARALLEL
		for (int i = 0; i < pairs.size(); i++)
			multiply(pairs[i].first, pairs[i].second);
	}




	template<typename group_t>
	void MultTable<group_t>::multiply_all_indices() {

		index_product.resize(degree.size(), basicIrreducibles.size());
		Greens.resize(degree.size());
		for (size_t i = 0; i < degree.size(); i++) {
			Greens[i].resize(basicIrreducibles.size());
			for (size_t j = 0; j < basicIrreducibles.size(); j++) {
				auto degreeproduct = degree[i] + basicIrreducibles[j];
				auto iterator = antidegree.find(degreeproduct);
				if (iterator == antidegree.end()) {
					index_product(i, j) = -1;
					continue;
				}
				index_product(i, j) = iterator->second;
			}
		}
	}

	///Compute the product of generators a_i*b_j
	template<typename group_t>
	void MultTable<group_t>::multiply(int i, int j) {
		auto C = basicChains[j];
		auto degreeC = Reindex<group_t>(basicIrreducibles[j]).front();
		auto degreeD = Reindex<group_t>(degree[i]).front();
		auto D = sphereChains.at(perfect_hasher(getsphere(degree[i])));
		Greens[i][j] = Green<group_t>(C, D, level, degreeC, degreeD);
		return;
	}

	template<typename group_t>
	Green<group_t> MultTable<group_t>::triple_product(int i, int j1, int j2) const {
		auto degreeC1 = Reindex<group_t>(basicIrreducibles[j1]).front();
		auto degreeC2 = Reindex<group_t>(basicIrreducibles[j2]).front();
		auto degreeD = Reindex<group_t>(degree[i]).front();
		auto C1 = basicChains[j1];
		auto C2 = basicChains[j2];
		auto D = ROChains<group_t>(getsphere(degree[i]));
		auto Box1 = Tensor<chains_t<group_t>>(C1, D, degreeC1 + degreeC2 + degreeD + 1).tensor();
		return Green<group_t>(C2, Box1, level, degreeC2, degreeC1 + degreeD);
	}
}
