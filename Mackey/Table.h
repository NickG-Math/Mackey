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

		///Extends multiplication table using the larger range provided
		void extend(const std::vector<int>&, const std::vector<int>&);

	protected:
		///////////////////////////////////////
		///The nonzero homology groups in our given range and level and their identification information (if non cyclic). 

		///The index i group corresponds to the antidegree[i] degree.
		///////////////////////////////////////
		std::vector<IDGenerators<rank_t>> NonZeroHomology;

		std::vector<std::vector<int>> degree; ///< Maps each index of NonZeroHomology to the corresponding degree
		std::map<std::vector<int>, int> antidegree; ///< Maps each degree to the corresponding index of NonZeroHomology

		Eigen::Matrix<int, -1, -1> index_product; ///<Maps each index and irreducible to the index of their product

		std::vector<int> minsphere;///<The lower bound on the range of our spheres
		std::vector<int> maxsphere;///<The upper bound on the range of our spheres

		std::vector<std::vector<int>> basicIrreducibles;///<The basic irreducibles we use to produce the factorizations (eg Euler and orientation classes)
		int number_of_irreducibles; ///<The number of basic irreducibles

		std::vector<std::vector<Green<rank_t, diff_t>>> Greens; ///<Each entry in the multiplication table
		std::map<std::array<int, 3>, Green<rank_t, diff_t> > tripleGreens; ///<The tripleGreens are set and used by triple_product

		/// Isolate sphere from the degree vector (i.e. get everything but first entry)
		inline std::vector<int> getsphere(const std::vector<int>& deg) {
			return tail(deg.data(), deg.size(),1);
		}

		/// Multiplies the given index with two basic irreducibles. Only used by MultiplicationGraphIdentify when all other identification techniques fail
		void triple_product(int, int, int);

#ifdef CEREALIZE
	public:
		/// Default constructor needed for loading archive
		MultiplicationTable() {};
		template<typename Archive, typename s_rank, typename s_diff>
		friend void serialize(Archive&, MultiplicationTable<s_rank, s_diff>&);
#endif

	private:
		std::vector<std::pair<int, int>> within_range;
		std::vector<Chains<rank_t, diff_t>> IndexedChains; ///<The Chains of the nonzero homology groups.
		std::vector<Chains<rank_t, diff_t>> basicChains;///<The Chains of the basic irreducibles.

		///Hashes the given sphere to produce the index for IndexedChains.
		template<typename sphere_t>
		int hash(const sphere_t& sphere) { return hashvector(sphere, minsphere, maxsphere); }

		void getIrreducibleChains();
		void multiply_all_indices();
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
	void MultiplicationTable<rank_t, diff_t>::getIrreducibleChains() {
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

		getIrreducibleChains();
		generators();
		multiply_all_indices();

#pragma omp parallel for num_threads(12) schedule(dynamic)
		for (int i = 0; i < within_range.size(); i++) {
			multiply(within_range[i].first, within_range[i].second);
		}
	}


	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::extend(const std::vector<int>& min, const std::vector<int>& max) {
		MultiplicationTable<rank_t, diff_t> M;
		M.level = level;
		M.minsphere = min;
		M.maxsphere = max;
		M.basicIrreducibles = std::move(basicIrreducibles);
		M.basicChains = std::move(basicChains);
		M.number_of_irreducibles = number_of_irreducibles;
		M.generators();
		M.multiply_all_indices();

		for (auto it = tripleGreens.begin(); it != tripleGreens.end(); it++)
			M.tripleGreens[{M.antidegree[degree[it->first[0]]], it->first[1], it->first[2]}] = std::move(it->second);

		std::vector<int> glossary;
		glossary.assign(M.degree.size(), -1);
		for (int i = 0; i < degree.size(); i++) {
			M.Greens[M.antidegree[degree[i]]] = std::move(Greens[i]);
			glossary[M.antidegree[degree[i]]] = i;
		}

		std::vector<std::pair<int, int>> to_be_done;
		to_be_done.reserve(M.degree.size() * number_of_irreducibles);
		for (int i = 0; i < M.degree.size(); i++) {
			for (int j = 0; j < number_of_irreducibles; j++) {
				if (M.index_product(i, j) != -1 && (glossary[i]==-1 || index_product(glossary[i], j) == -1) ) //if not original degree, or maybe out of range for original
					to_be_done.push_back(std::make_pair(i, j));
			}
		}

#pragma omp parallel for num_threads(12) schedule(dynamic)
		for (int i = 0; i < to_be_done.size(); i++) 
			M.multiply(to_be_done[i].first, to_be_done[i].second);

		*this = M;
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::multiply_all_indices() {

		index_product.resize(degree.size(), number_of_irreducibles);
		within_range.reserve(index_product.size());
		Greens.resize(degree.size());

		for (int i = 0; i < degree.size(); i++) {
			Greens[i].resize(number_of_irreducibles);
			for (int j = 0; j < number_of_irreducibles; j++) {
				auto degreeproduct = degree[i] + basicIrreducibles[j];
				auto iterator = antidegree.find(degreeproduct);
				if (iterator == antidegree.end()) {
					index_product(i, j) = -1;
					continue;
				}
				index_product(i, j) = iterator->second;
				within_range.push_back(std::make_pair(i, j));
			}
		}
	}

	///Compute the product of generators a_i*b_j
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::multiply(int i, int j) {
		auto C = basicChains[j];
		auto D = IndexedChains[hash(getsphere(degree[i]))];
		auto degreeC = Reindex(basicIrreducibles[j]).front();
		auto degreeD = Reindex(degree[i]).front();
		Greens[i][j] = Green<rank_t, diff_t>(C, D, level, degreeC, degreeD);
		return;
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::triple_product(int i, int j1, int j2) {
		auto i_j1 = index_product(i, j1); //should be in range
		auto i_j1_j2 = index_product(i_j1, j2); //should be in range
		auto degreeC1 = Reindex(basicIrreducibles[j1]).front();
		auto degreeC2 = Reindex(basicIrreducibles[j2]).front();
		auto degreeD = Reindex(degree[i]).front();
		auto C1 = basicChains[j1];
		auto C2 = basicChains[j2];
		auto D = IndexedChains[hash(getsphere(degree[i]))];
		auto Box1 = Box(C1, D, degreeC1 + degreeC2 + degreeD + 1);
		tripleGreens[{i, j1, j2}] = Green<rank_t, diff_t>(C2, Box1, level, degreeC2, degreeC1 + degreeD);
		return;
	}
}
