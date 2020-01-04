#pragma once
#include "Chains.h"
#include "General.h"
#include "Smith.h"
#include "Z_n.h"
#include <numeric>

///@file
///@brief Contains the Homology class and algorithms.

namespace Mackey {
	

	///The Homology of a Junction
	template<typename rank_t, typename diff_t>
	class Homology {

		///<Type used for Smith + product of matrices, either for performance (float/double) or to avoid overflow (long).
		typedef typename std::conditional<is_Finite_Cyclic<Scalar_t<diff_t>>::value, Scalar_t<diff_t>, float>::type fScalar;
		typedef typename std::conditional<is_Dense<diff_t>::value, Eigen::Matrix<fScalar, -1, -1>, Eigen::SparseMatrix<fScalar, 0>>::type fdiff_t_C;	///<Column Major
		typedef typename std::conditional<is_Dense<diff_t>::value, Eigen::Matrix<fScalar, -1, -1, 1>, Eigen::SparseMatrix<fScalar, 1>>::type fdiff_t_R;	///<Row Major
		typedef typename std::conditional<is_Dense<diff_t>::value, Eigen::Matrix<fScalar, -1, 1>, Eigen::SparseVector<fScalar>>::type fdiff_t_C_col;	///<Column of column major
	
	public:

		///The type of our matrix of generators
		typedef fdiff_t_C Gen_t;
		///The type of our generators (a column in the generator matrix)
		typedef fdiff_t_C_col gen_t;


		rank_t Groups;///<Encodes the homology groups as follows: Groups=[1,2,3] means homology Z+Z/2+Z/3
		Gen_t Generators;///<Encodes the generators homology groups as follows: The i-th column corresponds to the generator for Groups[i]
		bool isZero;///<1 if the homology is trivial
		
		///Default constructor
		Homology() {};

		///Compute the homology from the given Junction
		Homology(const Junction<rank_t, diff_t>&);

		///Compute the homology from the given Junction with an extra computation of a coefficient matrix needed for boundaries (if getQ=1)
		Homology(const Junction<rank_t, diff_t>&, bool);

		//////////////////////////////////////
		///Given an element in homology, write it as a linear combination of the generators of the homology

		///The answer is encoded as follows: basis=[-1,0,3] means element=-gen[0]+3*gen[2]
		///////////////////////////////////////
		rank_t basis(const gen_t&) const;

		///Given an x that is a boundary returns a y s.t. dy=x
		gen_t boundary(const gen_t&) const;

	public:
		std::vector<int> dontModOut;
		fdiff_t_C In_Q;
		fdiff_t_R Out_Qi, In_P_full, In_P_reduced;
		row_t<fdiff_t_C> diagonal;
		int M;
		fdiff_t_C getKernel(fdiff_t_C&);
		void KernelModImage(fdiff_t_C&, fdiff_t_C&, bool);

#ifdef CEREALIZE
		template<typename Archive, typename srank_t, typename sdiff_t>
		friend void serialize(Archive&, Homology<srank_t, sdiff_t>&);
#endif
	};


	template<typename rank_t, typename diff_t>
	Homology<rank_t, diff_t>::Homology(const Junction<rank_t, diff_t>& J, bool getQ) {
		M = summation(J.rank);
		fdiff_t_C In, Out;
		if (J.diffIn.size() == 0) {
			In.resize(M, 1);
			In.setZero();
		}
		else
			In= J.diffIn.template cast<fScalar>();
		if (J.diffOut.size() == 0) {
			Out.resize(1,M);
			Out.setZero();
		}
		else 
			Out = J.diffOut.template cast<fScalar>();
		auto Kernel=getKernel(Out);
		if (isZero)
			return;
		KernelModImage(In, Kernel, getQ);
	}

	template<typename rank_t, typename diff_t>
	Homology<rank_t, diff_t>::Homology(const Junction<rank_t, diff_t>& J) :Homology(J, 0) {}

	template<typename rank_t, typename diff_t>
	typename Homology<rank_t, diff_t>::fdiff_t_C Homology<rank_t, diff_t>::getKernel(fdiff_t_C& Out) {
		auto OUT = diagonalize<fdiff_t_C, fdiff_t_R, fdiff_t_C>(Out, 0, 1);
		fdiff_t_C Kernel(M,M);
		std::vector<int> nonZeroVectors;
		nonZeroVectors.reserve(M);
		isZero = 1;
		int j = 0;
		for (int i = 0; i < M;i++) {
			if (i >= OUT.diagonal.size() || OUT.diagonal[i] == 0) {
				isZero = 0;
				Kernel.col(j) = OUT.Q.col(i);
				nonZeroVectors.push_back(i);
				j++;
			}
		}
		Out_Qi = KeepRow(OUT.Qi, nonZeroVectors);
		if (!isZero)
			Kernel.conservativeResize(M, j);
		return Kernel;
	}

	template<typename rank_t, typename diff_t>
	void Homology<rank_t, diff_t>::KernelModImage(fdiff_t_C& In, fdiff_t_C& Kernel, bool getQ) {
		In = Out_Qi * In;
		if constexpr (is_Sparse<diff_t>::value)
			In.pruned();
		auto L = std::min(In.rows(), In.cols());
		auto IN=diagonalize<fdiff_t_C, fdiff_t_R, fdiff_t_C>(In, 1, getQ, 1);
		Generators = Kernel * IN.Pi;
		if constexpr (is_Sparse<diff_t>::value)
			Generators.pruned();
		auto maxsize = Generators.cols();
		std::vector<Scalar_t<rank_t>> groups;
		groups.reserve(maxsize);
		dontModOut.reserve(maxsize);
		isZero = 1;
		diagonal = std::move(IN.diagonal);
		for (int i = 0; i < maxsize;i++) {
			if (i < L) {
				if (abs(diagonal[i]) == 0) {
					groups.push_back(1);
					isZero = 0;
					dontModOut.push_back(i);
				}
				else if (abs(diagonal[i]) != 1) {
					groups.push_back(static_cast<Scalar_t<rank_t>>(abs(diagonal[i])));
					isZero = 0;
					dontModOut.push_back(i);
				}
			}
			else {
				groups.push_back(1);
				isZero = 0;
				dontModOut.push_back(i);
			}
		}
		In_P_reduced = KeepRow(IN.P,dontModOut);
		if (getQ) {
			In_P_full = std::move(IN.P);
			In_Q = std::move(IN.Q);
		}

		if (isZero) {
			Generators.resize(0,0);
			return;
		}
		Groups = Eigen::Map<rank_t>(groups.data(), groups.size());
		Generators = KeepCol(Generators, dontModOut);


		//Check for non integer coefficients and replace the Z in Groups with Z/N
		if constexpr (is_Finite_Cyclic<Scalar_t<diff_t>>::value) {
			constexpr int order = Scalar_t<diff_t>::order;
			for (int i = 0; i < Groups.size(); i++) {
				if (Groups[i] == 1)
					Groups[i] = static_cast<Scalar_t<rank_t>>(order);
			}
		}
	}


	template<typename rank_t, typename diff_t>
	rank_t Homology<rank_t, diff_t>::basis(const gen_t& generator) const {
		if (isZero)
			return rank_t();
		rank_t basisArray(Groups.size());
		gen_t element = In_P_reduced * Out_Qi * generator;
		for (int j = 0; j < Groups.size();j++) {
			if (Groups[j] != 1)
				basisArray[j] = static_cast<Scalar_t<rank_t>>((Groups[j] + (long)element[j] % Groups[j]) % Groups[j]);
			else
				basisArray[j] = static_cast<Scalar_t<rank_t>>(element[j]);
		}
		return basisArray;
	}



	template<typename rank_t, typename diff_t>
	typename Homology<rank_t, diff_t>::gen_t Homology<rank_t, diff_t>::boundary(const gen_t& generator) const {
		gen_t element = In_P_full * Out_Qi * generator;
		for (const auto& i : dontModOut) {
			if (element[i] != 0)  //the element is not 0 so it has no preimage
				return gen_t();
		}
		gen_t y(In_Q.rows()); //Sy=Px
		y.setZero();
		for (int i = 0; i < y.size(); i++) {
			if (i < diagonal.size() && diagonal[i] != 0)
				y[i] = static_cast<fScalar>(element[i] / diagonal[i]);
		}
		return In_Q * y;
	}


	////////////////////////////////////////////////////
///Finds the order of element in given finitely generated abelian group. Returns 0 if order is infinite

///T,S are vectors or 1d Eigen matrices
///////////////////////////////////////////////////
	template<typename T, typename S>
	int order(const T& element, const S& group) {
		std::vector<int> fractions;
		fractions.reserve(group.size());

		for (int i = 0; i < group.size(); i++) {
			if (group[i] == 1 && element[i] != 0)
				return 0; //infinite order
			else if (element[i] != 0) {
				fractions.push_back(group[i] / std::gcd(element[i], group[i]));
			}
		}
		return lcm(fractions);
	}

	///Normalizes the element in a group so that 0<=element[i]<=group[i] if group[i]!=1
	template<typename T>
	void mod_normalize(T& element, const T& group) {
		for (int i = 0; i < group.size(); i++) {
			if (group[i] != 1)
				element[i] = (element[i] + group[i]) % group[i];
		}
	}


	///Normalizes the element in a group to have minimal signs amongst its multiples
	template<typename T>
	void normalize(Eigen::Matrix<T,1,-1>& element, const Eigen::Matrix<T, 1, -1>& group) {
		if (element.isZero())
			return;
		mod_normalize(element,group);
		if (element.isZero())
			return;
		int size = element.size();
		if (size == 1) {
			if (group[0] == 1 && element[0] < 0) {
				element = -element;
				mod_normalize(element, group);
				return;
			}
			if (group[0]!=1) {
				element[0] = std::gcd(element[0], group[0]);
				return;
			}
		}
		auto n = order(element, group);
		if (n == 0) {
			for (int i = 0; i < size; i++) {
				if (group[i] == 1 && element[i] < 0) {
					element = -element;
					mod_normalize(element, group);
					return;
				}
			}
		}
		if (n <= 2)
			return;
		else {
			auto ideal_element = element;
			for (int i = 0; i < size; i++) {
				if (group[i] != 1 && ideal_element[i] != 0)
					ideal_element[i] = std::gcd(ideal_element[i], group[i]);
			}
			int closest_to_ideal = 1;
			int max_hits = -1;
			for (int i = 1; i < n; i++) {
				if (std::gcd(i, n) == 1) {
					int current_hits = 0;
					for (int j = 0; j < size; j++) {
						if ((i * element[j] - ideal_element[j]) % group[j]==0)
							current_hits++;
					}
					if (max_hits == -1 || max_hits < current_hits) {
						if (current_hits == size) {
							element = ideal_element;
							return;
						}
						max_hits = current_hits;
						closest_to_ideal = i;
					}
				}
			}
			if (closest_to_ideal == 1)
				return;
			element=closest_to_ideal * element;
			mod_normalize(element, group);
			return;
		}
	}


	////////////////////////////////////////////////////
///Finds if two elements are equal in the given finitely generated abelian group

///T,S,R are vectors or 1d Eigen matrices
///////////////////////////////////////////////////
	template<typename T, typename S, typename R>
	bool equals(const T& element1, const S& element2, const R& group) {
		for (int i = 0; i < group.size(); i++) {
			if (group[i] == 1 && element1[i] != element2[i])
				return 0;
			if (group[i] != 1 && (((element1[i] - element2[i]) % group[i]) != 0))
				return 0;
		}
		return 1;
	}

	////////////////////////////////////////////////////
///Finds if element a is a multiple of element b in the given finitely generated abelian group

///T,S,R are vectors or 1d Eigen matrices
///////////////////////////////////////////////////
	template<typename T, typename S>
	int isMultiple(const T& a, const T& b, const S& group) {
		auto oa = order(a, group);
		auto ob = order(b, group);
		if ((oa == 0 && ob != 0) || (ob < oa))
			return 0;
		if (oa == 0) {
			Scalar_t<T> q = 0;
			for (int i = 0; i < a.size(); i++) {
				if (group[i] == 1 && a[i] != 0) {
					if (b[i] == 0)
						return 0;
					else {
						q = a[i] / b[i];
						break;
					}
				}
			}
			if (equals(a, q * b, group))
				return q;
			else
				return 0;
		}
		for (Scalar_t<T> i = 1; i < ob; i++) {
			if (equals(a, i * b, group))
				return i;
		}
		return 0;
	}
	

	// //temp sparse implementation
	// template<typename rank_t, typename diff_t>
	// Homology<rank_t, diff_t>::Homology(const Junction<rank_t, diff_t>& J){
		// Junction<rank_t, spm_t<diff_t>> Jnew(J.rank, J.rankOut, J.rankIn, J.diffOut.sparseView(), J.diffIn.sparseView());
		// Homology<rank_t, spm_t<diff_t>> H(Jnew,0);
		// Groups = H.Groups;
		// isZero = H.isZero;
		// Generators = H.Generators;
		// Out_Qi = H.Out_Qi;
		// In_P_reduced = H.In_P_reduced;
	// };




// untested
//	////////////////////////////////////////////////////
/////Finds if element a is a linear combination of other given elements in a finitely generated abelian group
//
/////T,S are vectors or 1d Eigen matrices
/////////////////////////////////////////////////////
//	template<typename T, typename S>
//	bool inSpan(const T& a, const std::vector<T>& b, const S& group) {
//		std::vector<int> minbasis(b.size());
//		std::vector<int> maxbasis, tobegcded;
//		maxbasis.reserve(b.size());
//		for (int i=0; i<group.size(); i++){
//			if (group[i]==0){
//				for (const auto & j: b)
//					tobegcded.push_back(j[i]);
//				maxbasis.push_back(gcd(tobegcded));
//			}
//			else
//				maxbasis.push_back(group[i]-1);
//		}
//		auto U=DegreeConstruction(minbasis,maxbasis);
//		for (const auto & i:U){
//			T element=i[0]*b[0];
//			for (int j=1; j<b.size(); j++)
//				element+=i[j]*b[j];
//			normalize(element,group);
//			if (a==element)
//				return 1;
//		}
//		return 0;
//	}


}
