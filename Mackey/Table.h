#pragma once
#include "Compute.h"
#include "Green.h"

///@file
///@brief Contains the classes and methods to form the multiplication table/graph.

namespace Mackey {

	/// The input for the multiplication table
	template<typename rank_t, typename diff_t>
	class TableInput {
	protected:

		///////////////////////////////////////////////////
		/// The nonzero homology groups used for the table. 

		/// The order is not important (we use antimap to recover the degrees and generators).
		///////////////////////////////////////////////////
		std::vector<rank_t> NonZeroHomology;

		///////////////////////////////////////////////////
		/// The Chains of the nonzero homology groups.

		/// They are saved for performance.
		///////////////////////////////////////////////////
		std::vector<Chains<rank_t, diff_t>> IndexedChains;

		///////////////////////////////////////////////////
		/// antimap[i] reveals the generator of NonZeroHomology[i] 

		/// antimap[i] consists of the degree, the position of the generator (0 unless noncyclic) and a multiple (if we don't have a generator, but rather the multiple of one)
		/// Example: {5,4,-2,0,2} means that we have 2 times the 0-th generator in degree k=5, n=4, m=-2.
		///////////////////////////////////////////////////
		std::vector<std::vector<int>> antimap;

		std::map<std::vector<int>, int> map;///<The inverse of antimap

		const std::vector<int> maxsphere;///<The upper bound on the range of our spheres
		const std::vector<int> minsphere;///<The lower bound on the range of our spheres

		/// Constructs the input of the multiplication table. Currently only for C4 (and top level) but is easily extendible
		TableInput(const std::vector<int>& minsphere, const std::vector<int>& maxsphere);

		/// Checks if a given sphere is within the range of minsphere and maxsphere
		template<typename deg_t>
		bool withinrange(const deg_t&);

		///Hashes the given sphere to produce the index for IndexedChains.
		template<typename deg_t>
		int hash(const deg_t&);
	};


	template<typename rank_t, typename diff_t>
	TableInput<rank_t, diff_t>::TableInput(const std::vector<int>& minsphere, const std::vector<int>& maxsphere) :
		minsphere(minsphere), maxsphere(maxsphere)
	{
		auto spheres = DegreeConstruction(minsphere, maxsphere);
		auto maxlength = spheres.size() * (dimension(maxsphere - minsphere) + 1);
		NonZeroHomology.reserve(maxlength);
		antimap.reserve(maxlength);
		IndexedChains.resize(spheres.size());

		for (const auto& i : spheres) {
			auto C = ROChains<rank_t, diff_t>(i);
			for (int k = 0; k <= C.maxindex;k++) {
				Junction<rank_t, diff_t> J(C, k);
				Junction<rank_t, diff_t> J_top = transfer(J, power);
				Homology<rank_t, diff_t> H(J_top);
				if (!H.isZero) {
					for (int j = 0; j < H.Groups.size(); j++) {
						NonZeroHomology.push_back(H.Groups); // we want NonZeroHomology to be in sync with antimap and map, this won't happen too many times

						std::vector<int> generator;
						generator.reserve(3 + i.size());
						generator.push_back(invReindex(k, i));
						generator.insert(generator.end(), i.begin(), i.end());
						generator.push_back(j);
						generator.push_back(1);
						antimap.push_back(generator);
						map[generator] = antimap.size() - 1;
					}
				}
			}
			IndexedChains[hash(i)] = std::move(C);
		}
	}


	template<typename rank_t, typename diff_t>
	template<typename deg_t>
	bool TableInput<rank_t, diff_t>::withinrange(const deg_t& sphere) {
		for (size_t i = 0; i < minsphere.size(); i++) {
			if (sphere[i + 1] < minsphere[i] || maxsphere[i] < sphere[i + 1]) {
				return 0;
			}
		}
		return 1;
	}


	template<typename rank_t, typename diff_t>
	template<typename deg_t>
	int TableInput<rank_t, diff_t>::hash(const deg_t& sphere) {
		int hash = sphere[0] - minsphere[0];
		for (size_t i = 1; i < sphere.size(); i++) {
			hash += (sphere[i] - minsphere[i]) * (maxsphere[i] - minsphere[i] + 1);
		}
		return hash;
	}
}

namespace {
	using namespace Mackey;

	template<typename rank_t, typename diff_t>
	std::vector<Chains<rank_t, diff_t>> getChains(const std::vector<std::vector<int>>& basicIrreducibles) {
		std::vector<Chains <rank_t, diff_t>> basicChains;
		basicChains.reserve(basicIrreducibles.size());
		for (const auto& i : basicIrreducibles) {
			basicChains.push_back(ROChains<rank_t, diff_t>(remove0th(i.data(),i.size())));
		}
		return basicChains;
	}
}

namespace Mackey{




	/// The Multiplication Table
	template<typename rank_t, typename diff_t>
	class MultiplicationTable : public TableInput<rank_t, diff_t> {
	protected:
		int number_of_nodes;///<The total number of elements in the multiplication graph (including multiples of generators)
		int second_pass_nodes;///<The total number of generators in the multiplication graph

		///////////////////////////////////////////////
		///The edges of the multiplication graph

		///We connect a to ab for basic irreducible b, and ab to a if multiplication by b is an injection
		///////////////////////////////////////////////
		std::vector<std::vector<int>> edges;

		///////////////////////////////////////////////
		///The colors of the edges of the multiplication graph

		///The edge from a to ab is red and the one from ab to a (if it exists) is blue
		///////////////////////////////////////////////
		std::vector<std::vector<char>> colors;

		const std::vector<std::vector<int>> basicIrreducibles;///<The basic irreducibles we use to produce the factorizations (eg Euler and orientation classes)
		const std::vector<Chains<rank_t, diff_t>> basicChains;///<The Chains of the basic irreducibles.

		///Constructs the multiplication table given the maximum and minimum spheres and the basic irreducibles.
		MultiplicationTable(const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&);
	private:
		std::pair<rank_t, rank_t> makeproduct(int, int);
		void makeEdgeCyclic(const rank_t&, const rank_t&, int, int, int, int);
		void firstpass();
		void secondpass();


		std::pair<int, int> idgenerator(const rank_t&, const rank_t&);
		rank_t adjustHomology(const rank_t&, const rank_t&);
		bool injection(int, int);
		bool checkAndAssign(std::vector<int>&, int);


	};

	template<typename rank_t, typename diff_t>
	MultiplicationTable<rank_t, diff_t>::MultiplicationTable(const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles)
		:TableInput<rank_t, diff_t>(minsphere, maxsphere), basicIrreducibles(basicIrreducibles), basicChains(getChains<rank_t,diff_t>(basicIrreducibles)){
		number_of_nodes = TableInput<rank_t, diff_t>::antimap.size();
		edges.reserve((prime + 1) * number_of_nodes); //when we need to add extra nodes
		edges.resize(number_of_nodes);
		colors.reserve((prime + 1) * number_of_nodes);
		colors.resize(number_of_nodes);
		firstpass();
		secondpass();
	};

	///////////////////////////////////////////
	///Identify generators for noncyclic groups

	///First is n if it's uniquely specified as n times a generator
	///First is -n if it's n times a dangerous generator
	///Second is the position of the unique generator or a possible position of a dangerous generator
	///Right now only implemented for Z/2+Z.
	template<typename rank_t, typename diff_t>
	std::pair<int, int> MultiplicationTable<rank_t, diff_t>::idgenerator(const rank_t& basis, const rank_t& Groups) {

		if (Groups.size() == 2 && Groups[0] == 2 && Groups[1] == 1) {
			if (basis[0] == 1 && basis[1] == 0) { //torsion
				return std::make_pair(1, 0);
			}
			else if (basis[1] > 0) { //non-torsion
				return std::make_pair(-basis(1), 1);
			}
		}
		return std::make_pair(0, 0);
	}

	///Adjust the group a multiple of a generator lives in
	///Eg if we have a Z/4 then twice the generator lives in a Z/2=2(Z/4).
	template<typename rank_t, typename diff_t>
	rank_t MultiplicationTable<rank_t, diff_t>::adjustHomology(const rank_t& basis, const rank_t& Groups) {
		if (Groups.size() == 1) {
			if (Groups[0] == 1) {
				return Groups;
			}
			else {
				return Groups / basis[0];
			}
		}
	}

	///Checks if Z/a->Z/b is an injection (note: if a=1 then we mean Z as opposed to Z/1)
	template<typename rank_t, typename diff_t>
	bool MultiplicationTable<rank_t, diff_t>::injection(int a, int b) {
		return ((a == 1 && b == 1) || (a != 1 && a <= b));
	}

	///Checks if k exists in v and pushes back if not
	template<typename rank_t, typename diff_t>
	bool MultiplicationTable<rank_t, diff_t>::checkAndAssign(std::vector<int>& v, int k)
	{
		for (auto i : v) {
			if (i == k)
				return 0;
		}
		v.push_back(k);
		return 1;
	}

	///Connect a_i to a_i*b_j and a_i*b_j to a_i if we have an injection
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::makeEdgeCyclic(const rank_t& normalBasis, const rank_t& Groups, int i, int j, int position_range, int multiple) {

		auto position_domain = this->antimap[i][this->antimap[i].size() - 2];//the second to last is position
		auto degreesum = this->antimap[i] + basicIrreducibles[j];

		//extend degreesum to include gen information
		degreesum.push_back(position_range); //gen placement
		degreesum.push_back(multiple); // multiple

		if (multiple == 1) {
			//get integer hash
			auto k = this->map.at(degreesum);
			if (checkAndAssign(edges[i], k)) {// draw edge from i to k with red color
				colors[i].push_back(0);
			}

			//if we have an injection
			if (injection(this->NonZeroHomology[i](position_domain), this->NonZeroHomology[k](position_range))) {
				if (checkAndAssign(edges[k], i)) { //draw the inverse edge from k to i with blue color
					colors[k].push_back(1);
				}
			}
		}
		else if (multiple != 0) {
			auto iterator = this->map.find(degreesum);
			if (iterator == this->map.end()) {

				//new node for multiple generator

				auto k = this->antimap.size();
				this->antimap.push_back(degreesum);
				this->NonZeroHomology.push_back(adjustHomology(normalBasis, Groups));
				this->map[degreesum] = k;
				edges.resize(k + 1);
				edges[k].reserve(2 * basicIrreducibles.size());
				colors.resize(k + 1);
				colors[k].reserve(2 * basicIrreducibles.size());
				edges[i].push_back(k);
				colors[i].push_back(0);
				if (injection(this->NonZeroHomology[i](position_domain), this->NonZeroHomology[k](position_range)))
				{
					edges[k].push_back(i);
					colors[k].push_back(1);
				}
			}
			else {
				auto k = iterator->second;
				if (checkAndAssign(edges[i], k)) {
					colors[i].push_back(0);
				}
				if (injection(this->NonZeroHomology[i](position_domain), this->NonZeroHomology[k](position_range)))
				{
					if (checkAndAssign(edges[k], i)) {
						colors[k].push_back(1);
					}

				}
			}
		}
	}

	///////////////////////////////////////////////////////////////
	///Forms the first batch of edges of the multiplication table.

	///Only the generators (not their multiples) are treated in this first pass.
	///////////////////////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::firstpass() {
		second_pass_nodes = number_of_nodes;

		std::vector<std::pair<rank_t, rank_t>> thepair(number_of_nodes * basicIrreducibles.size());

//#pragma omp parallel for num_threads(12)
		for (int i = 0; i < number_of_nodes; i++) {
			edges[i].reserve(2 * basicIrreducibles.size());
			for (int j = 0; j < basicIrreducibles.size(); j++) {
				thepair[j + i * basicIrreducibles.size()] = makeproduct(i, j);
			}
		}
		for (int i = 0; i < number_of_nodes; i++) {
			for (int j = 0; j < basicIrreducibles.size(); j++) {
				auto normalBasis = thepair[j + i * basicIrreducibles.size()].first;
				auto Groups = thepair[j + i * basicIrreducibles.size()].second;
				if (normalBasis.size() == 0) {
					continue;
				}
				if (normalBasis.size() == 1)
				{
					makeEdgeCyclic(normalBasis, Groups, i, j, 0, normalBasis[0]);
				}
				else //non cyclic
				{
					auto a = idgenerator(normalBasis, Groups);
					if (a.first > 0) {
						makeEdgeCyclic(normalBasis, Groups, i, j, a.second, a.first);
					}
					else if (a.first < 0)
					{
						//for now do nothing
					}
				}
			}
		}
		number_of_nodes = this->antimap.size();
	}

	///////////////////////////////////////////////////////////////
	///Forms the second and final batch of edges of the multiplication table
	///
	///If a is a generator connected to c then 2*a is now connected to 2*c.
	///////////////////////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::secondpass() {
		for (int i = second_pass_nodes; i < number_of_nodes; i++) {

			//find what node it is the multiple of 
			auto original_degree = this->antimap[i];
			auto multiple = original_degree.back();
			original_degree.back() = 1;
			auto original_i = this->map.at(original_degree);

			//connect through multiplication
			edges[original_i].push_back(i);
			colors[original_i].push_back(0);

		}
	}

	///Compute the product of generators
	template<typename rank_t, typename diff_t>
	std::pair<rank_t, rank_t> MultiplicationTable<rank_t, diff_t>::makeproduct(int i, int j) {

		auto degreesum = this->antimap[i] + basicIrreducibles[j];
		if (!this->withinrange(degreesum)) { //outside of range so don't bother
			rank_t normalBasis, Groups;
			return std::make_pair(normalBasis, Groups);
		}

		std::vector<int> index((this->antimap[i]).begin() + 1, (this->antimap[i]).end() - 2);
		auto hashed = this->hash(index);
		auto degreeD = Reindex(this->antimap[i]).front();
		auto D = this->IndexedChains[hashed];
		auto selectD = this->antimap[i][this->antimap[i].size()-2];
		auto degreeC = basicIrreducibles[j].front();
		auto C = basicChains[j];
		Green<rank_t, diff_t> G(C, D, power, degreeC, degreeD, 0, selectD);
		return std::make_pair(G.normalBasis, G.Groups);
	}
}
