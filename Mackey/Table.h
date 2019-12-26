#pragma once
#include "Compute.h"
#include "Green.h"
///@file
///@brief Contains the class and methods to form the multiplication table.

namespace Mackey
{
	/// The Multiplication Table
	template<typename rank_t, typename diff_t>
	class MultiplicationTable {
	public:
		int level; ///<The Mackey functor level we are working in.

		///Constructs the multiplication table given the maximum and minimum spheres and the basic irreducibles.
		MultiplicationTable(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&);

	protected:
		///////////////////////////////////////
		///The nonzero homology groups in our given range and level and their identification information (if non cyclic). 

		///The index i group corresponds to the antidegree[i] degree.
		///////////////////////////////////////
		std::vector<IDGenerators<rank_t>> NonZeroHomology;

		std::vector<std::vector<int>> degree; ///< Maps each index of NonZeroHomology to the corresponding degree
		std::map<std::vector<int>, int> antidegree; ///< Maps each degree to the corresponding index of NonZeroHomology

		std::vector<int> minsphere;///<The lower bound on the range of our spheres
		std::vector<int> maxsphere;///<The upper bound on the range of our spheres

		std::vector<std::vector<int>> basicIrreducibles;///<The basic irreducibles we use to produce the factorizations (eg Euler and orientation classes)
		int number_of_irreducibles; ///<The number of basic irreducibles

		std::vector<std::pair<Green<rank_t, diff_t>, int>> Greens; ///<The result of the multiplication table

		/// Isolate sphere from the degree vector
		inline std::vector<int> getsphere(const std::vector<int>& deg) {
			return tail(deg.data(), deg.size());
		}

#ifdef CEREALIZE
	public:
/// Default constructor needed for loading archive
		MultiplicationTable() {};
		template<typename Archive, typename s_rank, typename s_diff>
		friend void serialize(Archive&, MultiplicationTable<s_rank, s_diff>&);
#endif

	private:
		std::vector<Chains<rank_t, diff_t>> IndexedChains; ///<The Chains of the nonzero homology groups.
		std::vector<Chains<rank_t, diff_t>> basicChains;///<The Chains of the basic irreducibles.

		///Hashes the given sphere to produce the index for IndexedChains.
		template<typename sphere_t>
		int hash(const sphere_t& sphere) { return hashvector(sphere, minsphere, maxsphere); }

		void getChains();
		void multiply(int, int);
		void make();
		void generators();
	};

	template<typename rank_t, typename diff_t>
	MultiplicationTable<rank_t, diff_t>::MultiplicationTable(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles)
		:level(level), minsphere(minsphere), maxsphere(maxsphere), basicIrreducibles(basicIrreducibles) {
		number_of_irreducibles = basicIrreducibles.size();
		make();
	};

	///Set the basicChains given the basicIrreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::getChains() {
		basicChains.reserve(number_of_irreducibles);
		for (const auto& i : basicIrreducibles)
			basicChains.push_back(ROChains<rank_t, diff_t>(getsphere(i)));
	}

	///Set the basicChains given the basicIrreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::generators() {
		auto spheres = DegreeConstruction(minsphere, maxsphere);
		auto maxlength = spheres.size() * (dimension(maxsphere - minsphere) + 1);
		NonZeroHomology.reserve(maxlength);
		degree.reserve(maxlength);
		IndexedChains.resize(spheres.size());

		int degindex = 0;
		for (const auto& i : spheres) {
			auto C = ROChains<rank_t, diff_t>(i);
			for (int k = 0; k <= C.maxindex;k++) {
				Junction<rank_t, diff_t> J(C, k);
				internal::IDGeneratorCompute<rank_t, diff_t>ID(level, J);
				if (ID.ID.group.size() != 0) {
					NonZeroHomology.push_back(ID.ID);
					std::vector<int> deg;
					deg.reserve(1 + i.size());
					deg.push_back(invReindex(k, i));
					deg.insert(deg.end(), i.begin(), i.end());
					degree.push_back(deg);
					antidegree[deg] = degindex;
					degindex++;
				}
			}
			IndexedChains[hash(i)] = std::move(C);
		}
	}

	///Compute the products of all elements and irreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::make() {
		getChains();
		generators();
		Greens.resize(degree.size() * number_of_irreducibles);
#pragma omp parallel for num_threads(12) schedule(dynamic)
		for (int index = 0; index < degree.size(); index++) {
			for (int j = 0; j < number_of_irreducibles; j++)
				multiply(index, j);
		}
	}


	/////Compute the products of the specific elements and irreducibles provided. Not currently used
	//template<typename rank_t, typename diff_t>
	//void MultiplicationTable<rank_t, diff_t>::make(const std::vector<int>& element_indices, const std::vector<int>& irreducible_indices) {
	//	getChains();
	//	std::vector<int> degree_indices;
	//	degree_indices.reserve(element_indices.size());
	//	for (const auto& i : element_indices) {
	//		if (degree_indices.empty() || degree_indices.back() != i)
	//			degree_indices.push_back(i);
	//	}
	//	auto maxdegree_plus1 = degree_indices.back() + 1; //it's the size if consecutive
	//	auto maxirr_plus1 = irreducible_indices.back() + 1;
	//	Greens.resize(maxdegree_plus1 * maxirr_plus1);

	//	#pragma omp parallel for num_threads(12) schedule(dynamic)
	//	for (int  q = 0; q < degree_indices.size(); q++) { //for OpenMP parallelization we don't use a range based loop
	//		for (const auto& j : irreducible_indices) {
	//			Greens[j + degree_indices[q] * maxirr_plus1] = multiply(degree_indices[q], j);
	//		}
	//	}
	//}

	///Compute the product of generators a_i*b_j
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::multiply(int index, int j) {

		auto degreeproduct = degree[index] + basicIrreducibles[j]; //we are using the degree index i here
		auto iterator = antidegree.find(degreeproduct);
		if (iterator == antidegree.end()) {
			Green<rank_t, diff_t> G;
			Greens[j + index * number_of_irreducibles] = std::make_pair(G, -1);
			return;
		}
		auto index_product = iterator->second;
		auto degreeD = Reindex(degree[index]).front();
		std::vector<int> sphere(degree[index].begin() + 1, degree[index].end());
		auto D = IndexedChains[hash(sphere)];

		auto degreeC = Reindex(basicIrreducibles[j]).front();
		auto C = basicChains[j];

		Green<rank_t, diff_t> G(C, D, level, degreeC, degreeD);
		Greens[j + index * number_of_irreducibles] = std::make_pair(G, index_product);
		return;
	}

}
