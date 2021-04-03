#pragma once
#include "Mackey_Functors/Levels.hpp"
#include "Mackey_Functors/Identify.hpp"
#include "Chain_Complexes/Box.hpp"

///@file
///@brief Contains the class and methods for multiplying generators.

namespace mackey
{

	template <typename group_t>
	typename group_t::rank_t multiply(const chains_t<group_t> &C, const chains_t<group_t> &D, int level, int degreeC, int degreeD, int selectC = 0, int selectD = 0);

	//forward declaration for Clang
	namespace internal
	{
		template <typename>
		class GreenCompute;
	}

	/// The result of multiplying generators in a Green functor
	template <typename group_t>
	class Green
	{
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

	public:
		AbelianGroup<rank_t> group;	 ///<The homology group the product lives in.
		bool isZero;				 ///<1 if the homology group the product lives in is 0.
		IDGenerators<rank_t> boxID;	 ///<Identifies the generators of the homology group the product lives in.

		///Sets to zero
		Green();

		///Computes the product of generators given chains, degrees and level
		Green(const Chains<rank_t, diff_t> &, const Chains<rank_t, diff_t> &, int, int, int);

		///Returns the (normalized) product of the given two elements (not necessarily generators)
		rank_t getNormalBasis(const rank_t &, const rank_t&) const;

		///Returns the (normalized) product of the selected generator with another element
		rank_t getNormalBasis(int i, const rank_t &b) const;

		///Returns the (normalized) product of the given element with the selected generator
		rank_t getNormalBasis(const rank_t&a, int j) const;

		///Returns the (normalized) product of the selected generators
		rank_t getNormalBasis(int i, int j) const;

		bool operator==(const Green<group_t> &) const;

		std::vector<rank_t> basis; ///<For each selection of generators the product is a linear combination of generators in its degree
		int first_number_selections;  ///<The number of selections of generators for the first factor
		int second_number_selections; ///<The number of selections of generators for the second factor

	private:

		///Selects the corrent basis and normalBasis index given selections.
		int select(int, int) const;

		template <typename>
		friend class internal::GreenCompute;
	};

	namespace internal
	{

		///Computes the homology at given level, the generators and their restrictions
		template <typename group_t>
		class ChainsLevelGen
		{
			typedef typename group_t::rank_t rank_t;
			typedef typename group_t::diff_t diff_t;

			typename Homology<rank_t, diff_t>::Gens_t Gens;
			gen_t<rank_t, diff_t> gen, res_gen;
			rank_t rank_level;
			bool isZero;
			ChainsLevelGen(const Chains<rank_t, diff_t> &C, int level, int degree);
			ChainsLevelGen(const Chains<rank_t, diff_t> &C, int level, int degree, int selection);

			/// Computes the product of generators in a Green functor
			template <typename>
			friend class GreenCompute;

			///Computes Massey products and their indeterminacy
			template <typename>
			friend class MasseyCompute;
		};

		///Computes the product of generators at the bottom level
		template <typename rank_t, typename diff_t>
		class ProductGen
		{
			gen_t<rank_t, diff_t> product;
			int padright, padleft;
			const Chains<rank_t, diff_t> &C, D;
			int degreeC, degreeD;
			ProductGen(const Chains<rank_t, diff_t> &C, const Chains<rank_t, diff_t> &D, int degreeC, int degreeD);

			void pad(const std::vector<int64_t> &detailedrank);

			void multiply(gen_t<rank_t, diff_t> gen1, gen_t<rank_t, diff_t> gen2);

			///Computes products
			template <typename>
			friend class GreenCompute;

			///Computes Massey products
			template <typename>
			friend class MasseyCompute;
		};

		/// Computes the product of generators in a Green functor
		template <typename group_t>
		class GreenCompute
		{
			typedef typename group_t::rank_t rank_t;
			typedef typename group_t::diff_t diff_t;

			Green<group_t> G; ///<The answer of the computation
			gen_t<rank_t, diff_t> product;
			const Chains<rank_t, diff_t> &C, D;
			typename Homology<rank_t, diff_t>::Gens_t GenC, GenD;
			rank_t rank_bottom, rank_level;
			const int level, degreeC, degreeD;
			Homology<rank_t, diff_t> H_Level;
			ChainsLevelGen<group_t> C_Variables, D_Variables;
			ProductGen<rank_t, diff_t> Restricted;

			///Computes the product of generators
			GreenCompute(const Chains<rank_t, diff_t> &C, const Chains<rank_t, diff_t> &D, int level, int degreeC, int degreeD);

			/// Computes the generators from their degrees and the selections
			void getGens();
			/// Computes the tensor product of Chains C,D
			void box();
			/// Computes the product of generators genC, genD by combining the results of getGens and box
			void compute();

			/// The result of multiplying generators in a Green functor
			template <typename>
			friend class mackey::Green;
		};
	}
}
#include "impl/Green.ipp"
