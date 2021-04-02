#pragma once
#include "Green.hpp"

///@file
///@brief Contains the class and methods for Massey Products.

namespace mackey {

	///Stores Massey products and their indeterminacy
	template<typename group_t>
	class Massey {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

	public:
		rank_t basis; ///<Expresses the Massey product as a linear combination of the generators of the box product.
		rank_t normalBasis; ///<Same as basis, but normalized
		AbelianGroup<rank_t> Groups; ///<The homology group the product lives in.
		std::array<rank_t, 2> indeterminacy; ///< The indeterminacy of the Massey product, stored as two bases
		IDGenerators<rank_t> boxID; ///< Identification data for the generators of the Massey product
		bool exists; ///< If the Massey product exists
		bool noIndeterminacy; ///<If the Massey product has no indeterminacy
		bool isZero; ///<If the Massey product is 0.

		///Set 0 by default
		Massey();

		///Compute Massey product given level, generator degrees and selections
		Massey(const chains_t<group_t>&, const chains_t<group_t>&, const chains_t<group_t>&, int, int, int, int, int, int, int, bool);
	};


	namespace internal
	{
		///Computes Massey products and their indeterminacy
		template<typename group_t>
		class MasseyCompute {
			typedef typename group_t::rank_t rank_t;
			typedef typename group_t::diff_t diff_t;

			Massey<group_t> Mass;

			const chains_t<group_t>& C, D, E;
			gen_t<rank_t, diff_t> resgenC, resgenD, resgenE, boundaryCD, boundaryDE;

			chains_t<group_t> BoxCD, BoxDE;
			Junction<rank_t, diff_t> Box;

			IDGeneratorCompute<group_t> ID;

			std::vector<std::vector<int64_t>> CD_detailedrank, DE_detailedrank;
			std::vector<int64_t> CD_E_detailedrank, C_DE_detailedrank;

			const int level, degreeC, degreeD, degreeE, selectC, selectD, selectE;
			Eigen::PermutationMatrix<-1, -1, int64_t> BoxDE_to_BoxCD;

			typedef typename mackey::Tensor<chains_t<group_t>, std::vector<std::vector<int64_t>>> tensortype;
			gen_t<rank_t, diff_t> box_boundary(const chains_t<group_t>&, const chains_t<group_t>&, const tensortype&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);
			gen_t<rank_t, diff_t> product_bottom(const chains_t<group_t>&, const chains_t<group_t>&, const std::vector<int64_t>&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);

			MasseyCompute(const chains_t<group_t>&, const chains_t<group_t>&, const chains_t<group_t>&, int, int, int, int, int, int, int, bool);

			void getGens();
			void box();
			void compute();
			void indeterminacy();

			///Stores Massey products and their indeterminacy
			template<typename>
			friend class mackey::Massey;
		};
	}
}
#include "impl/Massey.ipp"