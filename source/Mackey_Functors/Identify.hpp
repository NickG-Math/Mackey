#pragma once
#include "Homology/Homology.hpp"
#include "Mackey_Functors/MackeyFunctor.hpp"

///@file
///@brief Contains the identification classes


namespace mackey {

	///Two levels of a Mackey functor, used for identification
	template<typename rank_t>
	class IDGenerators {
	public:
		///The matrix type of the transfer and restriction matrices
		typedef Eigen::Matrix<typename rank_t::Scalar, -1, -1> matrix_t;
		AbelianGroup<rank_t> group; ///<The homology group at our level
		AbelianGroup<rank_t> group_lower; ///<The homology group at one level lower
		matrix_t Tr; ///<The transfer from one level lower to our level
		matrix_t Res; ///<The restriction from our level to one lower

		///Default Constructor
		IDGenerators()=default;

		///Get a two-level Mackey functor out of this data
		MackeyFunctor<rank_t> getMackey() const;

		///Standard equality tests if all group, group_lower, Tr, Res are equal 
		bool operator==(const IDGenerators<rank_t>&) const;
	};

	///	Use the ID data to identify all possible candidates for an element
	template<typename rank_t>
	std::vector<rank_t> id_candidates(const rank_t&, const IDGenerators<rank_t>&, const IDGenerators<rank_t>&);


	///	Find if given elements can all be distinguished i.e. they have pairwise disjoint candidate sets
	template<typename rank_t>
	bool distinguish(const std::vector<rank_t>&, const IDGenerators<rank_t>&);


	///Contains classes and methods whose user interface is provided by the rest of the library
	namespace internal {
		///Computes the homology at given level and the identification data coming from the Mackey functor structure
		template<typename group_t>
		class IDGeneratorCompute {
			typedef typename group_t::rank_t rank_t;
			typedef typename group_t::diff_t diff_t;
		public:
			IDGenerators<rank_t> ID; ///<The identification data
			IDGeneratorCompute(int, const Junction<rank_t, diff_t>&, bool=0); ///<Computes the identification data given level, Junction at the bottom, and optionally if we want to store the Q matrix.
		private:
			rank_t rank_level;
			Homology<rank_t, diff_t> H_level;
			IDGeneratorCompute() = default;

			//Identifies the generators in case of non cyclic homology.
			void id(const Junction<rank_t, diff_t>&);

			///Compactifies the ID data
			template<typename>
			friend class IDGenerators;

			///Forms the input of the Multiplication Table
			template<typename, typename>
			friend class TableInput;

			///Computes products
			template<typename>
			friend class GreenCompute;

			///Computes Massey products
			template<typename>
			friend class MasseyCompute;

		};
	}
}
#include "impl/Identify.ipp"