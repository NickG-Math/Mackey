#pragma once
#include "Mult_Table.hpp"
#include "Utility/Graph.hpp"

///@file
///@brief Contains the multiplication graph.


namespace mackey {

	namespace implementation_details {
		struct compareElements;
	}

	/////////////////////////////////////////////////////////
	/// The Multiplication Graph created from the Multiplication Table

	/// A node is an element in a nonzero homology group, a linear combination of the generators. At initialization we only use generators, but as products of generators might not be generators,
	/// after that initial stage we add multiples and linear combinations into the mix.
	/// We have a new index for the elements, and use various vectors and ordered maps to keep track of the element index and the degree index. 
	/////////////////////////////////////////////////////////
	template<typename group_t>
	class MultGraph : public MultTable<group_t> {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

	public:
		int number_of_generators;///<The number of generators in the multiplication graph

		/// Retrieve the degree of the i-th generator
		std::vector<int> getdegree(int i) const;

		/// Retrieve the element the i-th generator corresponds to. 
		Eigen::Matrix<int, 1, -1> getelement(int i) const;

		/// Retrieve the element index of the given degree and element. Returns -1 if no such degree can be found
		int getelementindex(const std::vector<int>& deg, const rank_t& elmnt) const;

		typedef Graph<Neighborhood<size_t, bicolored_edge_with_id, std::vector>> graph_t;
		graph_t graph;

		///Constructs the multiplication graph given the maximum and minimum spheres and the basic irreducibles.
		template<typename ...Args>
		MultGraph(Args&&...);

	protected:


		std::vector<rank_t> element; ///<An element (linear combination of generators) in each NonZeroHomology group of the table.
		std::vector<int> tracker; ///< Maps element index to degree index

		/// Maps (degree_index, element)->element_index. The comparator first compares the degree_index and if equal then compares the elements first to last entry
		std::map<std::pair<int, rank_t>, int, implementation_details::compareElements> antielement;

		///Returns k if element[i]=k*element[j]. Returns 0 if element[a] is not a multiple of element[b]
		int isMultiple(int i, int j) const;

	private:
		void make();
		std::vector<std::pair<int, int>> unidentified; ///All pairs (i,j) for which we couldn't identify the product i*j
		std::vector<std::vector<int>> zeroproduct; ///<zeroproduct[i] consists of all j for which i*j leads to 0.

		void initialize();

		///Adds all edges starting and ending from/to given element
		void add_edges(int);
		///Connect each element to its multiple
		void pass_multiples();

		///Finds the index of the element given degree and presentation and adds it to the list if it doesn't exist there.
		int find_and_add_element(int, const rank_t&);

		///Determines how and if i can be connected to j*i
		auto determine_connection(const Green<group_t>&, int, int, int);

		///Connects i to j*i
		void connect(const rank_t&, int, int, int);

		std::vector<int> other_elements(int) const;

		///Checks if the surjective map from the i'th to the k'th generators is an injection
		bool injection(int, int) const;

		///Finds the order of the i'th generator
		int order(int i) const;

		template<typename>
		friend class MultIdentify;

		template<typename>
		friend class MultConnectivity;
	};
}
#include "impl/Mult_Graph.ipp"
