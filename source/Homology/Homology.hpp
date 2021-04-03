#pragma once
#include "Chain_Complexes/Chains.hpp"
#include "Utility/General.hpp"
#include "Smith.hpp"
#include "Abelian.hpp"

///@file
///@brief Contains the Homology class and algorithms.

namespace mackey {

	///The Homology of a Junction
	//
	///@tparam rank_t The rank type is a dense Eigen row vector of signed integer scalar eg Eigen::Matrix<int,1,-1>. 
	///It stores the ranks of the modules of the differentials 
	///Eg in the G equivariant case \f$[1,2,4]\f$ means \f$Z[G/G}\oplus Z[G/G']\oplus Z[G/G'']\f$ where \f$|G:G'|=2, |G/G''|=4\f$ (current implementation only for prime power cyclic groups)
	///@tparam diff_t The different
	template<typename rank_t, typename diff_t>
	class Homology {

		typedef std::conditional_t<SFINAE::is_finite_cyclic<scalar_t<diff_t>>::value, scalar_t<diff_t>, int64_t> HScalar; //scalar with extended accuracy for Z coefficients (no change for Z/N coefficients)
		typedef std::conditional_t<SFINAE::is_Dense<diff_t>::value, Eigen::Matrix<HScalar, -1, -1>, Eigen::SparseMatrix<HScalar, 0, storage_t<diff_t>>> diff_t_C; //diff_t column major with the scalar above
		typedef std::conditional_t<SFINAE::is_Dense<diff_t>::value, Eigen::Matrix<HScalar, -1, -1, 1>, Eigen::SparseMatrix<HScalar, 1, storage_t<diff_t>>> diff_t_R; //diff_t row major with the scalar above
	public:

		///The type of our matrix of generators
		typedef diff_t_C Gens_t;
		///The dense type of our generators (a column in the generator matrix, always dense for convenience)
		typedef col_vector_t<diff_t_C> gen_t;

		AbelianGroup<rank_t> group;///<Encodes the homology groups as follows: group=[1,2,3] means homology Z+Z/2+Z/3. This works even for Z/n coefficients: the free module (Z/n)^, is encoded as [n,n,...,n]
		Gens_t generators;///<Encodes the generators homology groups as follows: The i-th column corresponds to the generator for group[i]
		bool isZero;///<1 if the homology is trivial
		
		///Default constructor
		Homology()=default;

		///Compute the homology of Junction from the given Junction
		///@param J The given Junction 
		///@param getQ Whether we want to store the Q matrix; this is used by the boundary function
		Homology(const Junction<rank_t, diff_t>& J, bool getQ=0);

		//////////////////////////////////////
		///Given an element in homology, write it as a linear combination of the generators of the homology

		///The answer is encoded as follows: basis=[-1,0,3] means element=-gen[0]+3*gen[2]
		///////////////////////////////////////
		rank_t basis(const gen_t&) const;

		///Given an x that is a boundary returns a y s.t. dy=x
		gen_t boundary(const gen_t&) const;

	private:
		diff_t_R Out_Qi, In_P_full, In_P_reduced;
		std::vector<typename diff_t::StorageIndex> dontModOut;
		diff_t_C In_Q;
		row_vector_t<diff_t_C> diagonal;
		typename diff_t::StorageIndex M;
		diff_t_C getKernel(diff_t_C&);
		void KernelModImage(diff_t_C&, diff_t_C&, bool);

	};
}
#include "impl/Homology.ipp"
