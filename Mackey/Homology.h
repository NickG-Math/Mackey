#pragma once
#include "Chains.h"
#include "General.h"
#include "Smith.h"

///@file
///@brief Contains the Homology class and algorithms.

namespace Mackey {
	

	///The Homology of a Junction
	template<typename rank_t, typename diff_t>
	class Homology {
	public:
		
//////////////////////////////////////////////////
///The floating point type for Homology computations

///Should be float to maximize performance, or potentially double for groups of order not a prime power.
/////////////////////////////////////////////////
		typedef typename diff_t::Scalar fScalar;
		typedef Eigen::Matrix<fScalar, -1, -1> fdiff_t;	///<The type of our casted differentials





		rank_t Groups;///<Encodes the homology groups as follows: Groups=[1,2,3] means homology Z+Z/2+Z/3
		diff_t Generators;///<Encodes the generators homology groups as follows: The i-th column corresponds to the generator for Groups[i]
		bool isZero;///<1 if the homology is trivial
		
		///Default constructor
		Homology() {};

		///Compute the homology from the given Junction
		Homology(const Junction<rank_t, diff_t>&);

		//////////////////////////////////////
		///Given an element in homology, write it as a linear combination of the generators of the homology

		///The answer is encoded as follows: basis=[-1,0,3] means element=-gen[0]+3*gen[2]
		///////////////////////////////////////
		template<typename gen_t>
		rank_t basis(const gen_t&);
	private:
		std::vector<int> nonZeroVectors;
		std::vector<int> dontModOut;
		fdiff_t Q0i, P1;
	};


	template<typename rank_t, typename diff_t>
	Homology<rank_t, diff_t>::Homology(const Junction<rank_t, diff_t>& J) {
		auto M = summation(J.rank);

		fdiff_t Out, In; //these two will change and J might be needed again, so do copy


		if (J.diffIn.size() == 0 && J.diffOut.size() == 0) {
			In = Eigen::MatrixBase<fdiff_t>::Zero(M, 1);
			Out = Eigen::MatrixBase<fdiff_t>::Zero(1, M);
		}
		else if (J.diffOut.size() == 0) {
			Out = Eigen::MatrixBase<fdiff_t>::Zero(1, M);
			In = (J.diffIn).template cast<fScalar>();
		}
		else if (J.diffIn.size() == 0) {
			In = Eigen::MatrixBase<fdiff_t>::Zero(M, 1);
			Out = J.diffOut.template cast<fScalar>();
		}
		else {
			Out = J.diffOut.template cast<fScalar>();
			In = J.diffIn.template cast<fScalar>();
		}
		Smith<fdiff_t, fdiff_t, fdiff_t> OUT(Out, 0, 1);
		Q0i = std::move(OUT.Qi); //OUT won't be needed outside of this call so don't copy

		fdiff_t Kernel(M, M);
		nonZeroVectors.reserve(M);
		isZero = 1;
		int j = 0;
		for (int i = 0; i < M;i++) {
			if (((i > OUT.L - 1) || (OUT.S(i, i) == 0))) {
				//Remember that Q is invertible so it can't have zero columns
				isZero = 0;
				Kernel.col(j) = OUT.Q.col(i);
				nonZeroVectors.push_back(i);
				j++;
			}
		}
		if (isZero) {
			return;
		}
		Kernel.conservativeResize(M, j);
		In = Q0i * In;
		In = KeepRow(In, nonZeroVectors);

		Smith<fdiff_t, fdiff_t, fdiff_t> IN(In, 1, 0);
		P1 = std::move(IN.P); //IN won't be needed outside of this call so don't copy
		Kernel = Kernel * IN.Pi;
		Generators = Kernel.template cast<typename diff_t::Scalar>();
		auto maxsize = Kernel.cols();
		Groups.resize(maxsize);
		dontModOut.reserve(maxsize);

		isZero = 1;
		for (int i = 0; i < maxsize;i++) {
			if (i < In.rows() && i < In.cols()) {
				if (IN.S(i, i) == 0) {
					Groups(i) = 1;
					isZero = 0;
					dontModOut.push_back(i);

				}
				else if (abs(IN.S(i, i)) != 1) {
					Groups(i) = static_cast<typename rank_t::Scalar>(abs(IN.S(i, i)));
					isZero = 0;
					dontModOut.push_back(i);
				}
			}
			else {
				Groups(i) = 1;
				isZero = 0;
				dontModOut.push_back(i);
			}
		}
		if (isZero) {
			Groups.resize(0);
			return;
		}

		Generators = KeepCol(Generators, dontModOut);
		Groups = KeepCol(Groups, dontModOut);
	}

	

	template<typename rank_t, typename diff_t>
	template<typename gen_t>
	rank_t Homology<rank_t, diff_t>::basis(const gen_t& generator) {
		if (isZero) {
			rank_t basisArray;
			return basisArray;
		}
		rank_t basisArray(Groups.size());
		Eigen::Matrix<fScalar,-1,1> element = generator.template cast<fScalar>();
		element = Q0i * element;
		element = KeepRow(element, nonZeroVectors);
		element = P1 * element;
		element = KeepRow(element, dontModOut);

		for (int j = 0; j < Groups.size();j++) {
			if (Groups(j) != 1) {
				basisArray(j) = (Groups(j) + static_cast<typename rank_t::Scalar>(element(j)) % Groups(j)) % Groups(j);
				//This is needed due to C++ conventions for % take a symmetric range w.r.t. 0 instead of >=0. So we can't just do element(j)%Groups(j)
				//Also we can't use our own mod since this is mod of floats over int (although we could just define mod to just cast to int.
			}
			else {
				basisArray(j) = static_cast<typename rank_t::Scalar>(element(j));
			}
		}
		return basisArray;
	}
}
