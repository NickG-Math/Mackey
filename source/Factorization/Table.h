#pragma once
#include "Compute.h"
#include "omp.h" 


///@file
///@brief Contains the class and methods to form the multiplication table.

namespace Mackey
{
	/// The Multiplication Table
	template<typename rank_t, typename diff_t>
	class MultiplicationTable {
	public:
		int level; ///<The Mackey functor level we are working in.

		///Constructs the multiplication table given the maximum and minimum spheres and the basic irreducibles. Optionally do the computation in batches (given number of batches) and save the results of the i-th batch in binary named "batch" << i
		MultiplicationTable(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, int=0, bool=0);

		///Extends multiplication table using the larger range provided. Optionally do the computation in batches (given number of batches) and save the results of the i-th batch in binary named "batch" << i
		void extend(const std::vector<int>&, const std::vector<int>&, int = 0, bool = 0);

		/// Retrieve the degree index of the given degree
		int getdegreeindex(const std::vector<int>&) const;

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
		std::map<std::array<int, 3>, Green<rank_t, diff_t> > tripleGreens; ///<The tripleGreens set and used by MultiplicationGraphIdentify when all other identification techniques fail

		/// Isolate sphere from the degree vector (i.e. get everything but first entry)
		inline std::vector<int> getsphere(const std::vector<int>& deg) {
			return tail(deg.data(), deg.size(), 1);
		}

		/// Multiplies the given index with two basic irreducibles. Only used by MultiplicationGraphIdentify when all other identification techniques fail
		Green<rank_t, diff_t> triple_product(int, int, int);

#ifdef CEREALIZE
	public:
		/// Default constructor needed for loading archive
		MultiplicationTable() {};
		template<typename Archive, typename s_rank, typename s_diff>
		friend void save(Archive&, const MultiplicationTable<s_rank, s_diff>&);

		template<typename Archive, typename s_rank, typename s_diff>
		friend void load(Archive&, MultiplicationTable<s_rank, s_diff>&);
#endif

	private:
		int num_spheres;
		std::vector<Chains<rank_t, diff_t>> basicChains;///<The Chains of the basic irreducibles.
		std::map<int, Chains<rank_t, diff_t>> sphereChains;
		void getIrreducibleChains();
		void multiply_all_indices();
		void multiply(int, int);
		void make(int=0, bool=0);
		void generators();
		void compute(const std::vector<std::pair<int, int>>&, int=0, bool=0);



		int hash(const std::vector<int>& sphere) {
			return hashvector(sphere, minsphere, maxsphere);
		}

		int hash(int i) {
			return hash(getsphere(degree[i]));
		}


	};

	template<typename rank_t, typename diff_t>
	MultiplicationTable<rank_t, diff_t>::MultiplicationTable(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles, int number_of_teams, bool serialize_each_step)
		:level(level), minsphere(minsphere), maxsphere(maxsphere), basicIrreducibles(basicIrreducibles) {
		number_of_irreducibles = basicIrreducibles.size();
		make(number_of_teams, serialize_each_step);
	};

	template<typename rank_t, typename diff_t>
	inline int MultiplicationTable<rank_t, diff_t>::getdegreeindex(const std::vector<int>& degree) const {
		auto iterator = this->antidegree.find(degree);
		if (iterator == this->antidegree.end())
			return -1;
		else
			return iterator->second;
	}

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
		std::map<std::vector<int>, IDGenerators<rank_t>> Homologymap;
		num_spheres = spheres.size();
#if defined(_OPENMP)
		//Make sure to lock if you enable multithreading here!
		omp_lock_t lock;
		omp_init_lock(&lock);
#endif
//#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic)
		//compute additive in parallel and save to map
		for (int i = 0; i < num_spheres; i++) {
			auto C = ROChains<rank_t, diff_t>(spheres[i]);
			std::vector<int> deg;
			deg.reserve(1 + spheres[i].size());
			for (int k = 0; k <= C.maxindex;k++) {
				Junction<rank_t, diff_t> J(C, k);
				internal::IDGeneratorCompute<rank_t, diff_t>ID(level, J);
				if (ID.ID.group.size() != 0) {
					deg.clear();
					deg.push_back(invReindex(k, spheres[i]));
					for (const auto& j : spheres[i])
						deg.push_back(j);
#if defined(_OPENMP)
					omp_set_lock(&lock);
#endif
					Homologymap[deg] = ID.ID;
#if defined(_OPENMP)
					omp_unset_lock(&lock);
#endif
				}
			}
		}
#if defined(_OPENMP)
		omp_destroy_lock(&lock);
#endif
		//reinterpret map into vector
		auto maxlength = num_spheres * (dimension(maxsphere - minsphere) + 1);
		NonZeroHomology.reserve(maxlength);
		degree.reserve(maxlength);
		int degindex = 0;
		for (auto it = Homologymap.begin(); it != Homologymap.end(); it++) {
			auto deg = it->first;
			NonZeroHomology.push_back(it->second);
			degree.push_back(deg);
			antidegree[deg] = degindex;
			degindex++;
		}

	}

	///Compute the products of all elements and irreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::make(int number_of_teams, bool serialize_each_step) {

		getIrreducibleChains();
		generators();
		multiply_all_indices();

		std::vector<std::pair<int, int>> within_range;
		within_range.reserve(index_product.size());
		for (int j = 0; j < number_of_irreducibles; j++) {
			for (int i = 0; i < degree.size(); i++)
				if (index_product(i, j) != -1)
					within_range.push_back(std::make_pair(i, j));
		}
		compute(within_range, number_of_teams, serialize_each_step);
	}


	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::extend(const std::vector<int>& min, const std::vector<int>& max, int number_of_teams, bool serialize_each_step) {

		std::vector<std::pair<int, int>> to_be_done;
		{
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
			to_be_done.reserve(M.degree.size() * number_of_irreducibles);
			for (int i = 0; i < M.degree.size(); i++) {
				for (int j = 0; j < number_of_irreducibles; j++) {
					if (M.index_product(i, j) != -1 && (glossary[i] == -1 || index_product(glossary[i], j) == -1)) //if not original degree, or maybe out of range for original
						to_be_done.push_back(std::make_pair(i, j));
				}
			}

			*this = M;
		}
		compute(to_be_done, number_of_teams, serialize_each_step);
	}



	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::compute(const std::vector<std::pair<int, int>>& pairs, int number_of_teams, bool serialize_each_step) {

		if (number_of_teams == 0) { //all in one go
			for (const auto& i : pairs) {
				auto sphere = getsphere(degree[i.first]);
				auto j = hash(sphere);
				if (sphereChains.count(j)==0)
					sphereChains[j] = ROChains<rank_t, diff_t>(sphere);
			}
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic)
			for (int i = 0; i < pairs.size(); i++)
				multiply(pairs[i].first, pairs[i].second);
		}

		else { //do in batches
			std::vector<int> weights(num_spheres);
			std::vector<std::vector<int>> items(num_spheres);
			for (int i = 0; i < pairs.size(); i++) {
				auto sphere = getsphere(degree[pairs[i].first]);
				weights[hash(pairs[i].first)]++;
				items[hash(pairs[i].first)] = sphere;
			}
			auto teams = equidistribute(weights, items, number_of_teams);
			std::vector<std::vector<std::pair<int, int>>> batches;
			batches.resize(number_of_teams);
			for (int i = 0; i < pairs.size(); i++) {
				for (int j = 0; j < teams.size(); j++)
					if (find(teams[j], hash(pairs[i].first)) != -1) {
						batches[j].push_back(pairs[i]);
						break;
					}
			}
			for (int i = 0; i < batches.size(); i++) {
				sphereChains.clear();
				for (const auto& j : teams[i]) {
					if (!items[j].empty())
						sphereChains[hash(items[j])] = ROChains<rank_t, diff_t>(items[j]);
				}
#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic)
				for (int j = 0; j < batches[i].size(); j++)
					multiply(batches[i][j].first, batches[i][j].second);
#ifdef CEREALIZE
				if (serialize_each_step)
					saver(*this, "batch" + std::to_string(i), "binary");
#endif
			}
		}
	}




	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::multiply_all_indices() {

		index_product.resize(degree.size(), number_of_irreducibles);
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
			}
		}
	}

	///Compute the product of generators a_i*b_j
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::multiply(int i, int j) {
		auto C = basicChains[j];
		auto degreeC = Reindex(basicIrreducibles[j]).front();
		auto degreeD = Reindex(degree[i]).front();
		auto D = sphereChains.at(hash(i));
		Greens[i][j] = Green<rank_t, diff_t>(C, D, level, degreeC, degreeD);
		return;
	}

	template<typename rank_t, typename diff_t>
	Green<rank_t, diff_t> MultiplicationTable<rank_t, diff_t>::triple_product(int i, int j1, int j2) {
		auto degreeC1 = Reindex(basicIrreducibles[j1]).front();
		auto degreeC2 = Reindex(basicIrreducibles[j2]).front();
		auto degreeD = Reindex(degree[i]).front();
		auto C1 = basicChains[j1];
		auto C2 = basicChains[j2];
		auto D = ROChains<rank_t, diff_t>(getsphere(degree[i]));
		auto Box1 = Box(C1, D, degreeC1 + degreeC2 + degreeD + 1);
		return Green<rank_t, diff_t>(C2, Box1, level, degreeC2, degreeC1 + degreeD);
	}
}
