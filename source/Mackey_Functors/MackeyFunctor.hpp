#pragma once
#include <string>
#include <iostream>
#include <numeric>
#include "Homology/Abelian.hpp"

///@file
///@brief Contains the class describing Mackey Functors.

namespace mackey
{
	///Isomorphism type of Mackey Functors
	template<typename rank_t>
	using iso_t = std::vector<dense_t<rank_t>>;

	///A Mackey Functor
	template<typename rank_t>
	class MackeyFunctor {
	public:

		std::vector<AbelianGroup<rank_t>> levels;///<The levels starting at level 0 (see Homology for the labeling)

		//////////////////////////////////////////////////////////////////////////////
		///The transfers from level i to level i+1. They are encoded as matrices because the levels may be noncyclic

		///Example (C4): If tr[0]=[2] and tr[1]=[1,2;1,3] then Tr_0^2(gen)=2*gen and Tr_2^4(gen0)=gen0+gen1 and Tr_2^4(gen1)=2*gen0+3*gen1
		//////////////////////////////////////////////////////////////////////////////
		std::vector<dense_t<rank_t>> tr;

		std::vector<dense_t<rank_t>> res;///<Same as transfers, but now using restrictions
		std::vector<dense_t<rank_t>> act;///<Same as transfers, but now using the act group action

		std::string name; ///<The name of the Mackey functor

		///Prints a Mackey functor ignoring the name
		std::string print() const;

		///Resize all member variables.
		void resize(int);

		///Default constructor
		MackeyFunctor() = default;

		///When there is only one generator at each level we can use this convenient constructor	
		MackeyFunctor(const rank_t&, int);

		///Compares two Mackey and returns 1 if equal.
		bool operator==(const MackeyFunctor<rank_t>&) const;

		///Applies given isomorphism to Mackey functor (the inverse of the isomorphism must also be provided)
		MackeyFunctor<rank_t> apply(const iso_t<rank_t>&, const iso_t<rank_t>&) const;

		///Returns all automorphisms of a Mackey Functor and their inverses
		std::pair<std::vector<iso_t<rank_t>>, std::vector<iso_t<rank_t>>> automorphisms() const;

		///Returns the isomorphism class of a Mackey Functor
		std::vector<MackeyFunctor<rank_t>> isomorphism_class() const;

		///Returns 1 if the Mackey functor belongs to the provided isomorphism class
		bool isomorphic(const std::vector<MackeyFunctor<rank_t>>& isoclass) const;

		///Returns 1 if the Mackey functors are isomorphic.
		bool isomorphic(const MackeyFunctor<rank_t>&);

		////////////////////////////////////////////////////////////////////////////////////////////
		///Computes name for certain Mackey functors in our universal notation

		///Eg 124 means groups Z, Z/2, Z/4 bottom to top, with transfers multiplication by 2 and restrictions multiplication by 1.
		///Eg 112 #012 means groups Z, Z, Z/2 bottom to top, with transfers multiplication by 1 and restrictions multiplication by 2 (the sharp signifies restriction and transfer swap)
		/////////////////////////////////////////////////////////////////////////////////////////
		void notation();
	};

	///Print Mackey Functor to output stream
	template<typename rank_t>
	std::ostream& operator<<(std::ostream&, const MackeyFunctor<rank_t>&);
}
#include "impl/MackeyFunctor.ipp"
