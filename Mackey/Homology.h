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
	public:	
		rank_t Groups;///<Encodes the homology groups as follows: Groups=[1,2,3] means homology Z+Z/2+Z/3
		diff_t Generators;///<Encodes the generators homology groups as follows: The i-th column corresponds to the generator for Groups[i]
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
		template<typename gen_t>
		rank_t basis(const gen_t&) const;

		///Given an x that is a boundary returns a y s.t. dy=x
		template<typename gen_t>
		gen_t boundary(const gen_t&) const;

	private:
		typedef float fScalar;	///<The floating point type used when taking products of matrices (and for the Smith when using Z coefficients)
		typedef Eigen::Matrix<fScalar, -1, -1> fdiff_t;	///<The type of our casted differentials

		std::vector<int> nonZeroVectors;
		std::vector<int> dontModOut;
		fdiff_t Out_Qi, In_P, In_Q;
		Eigen::Matrix<typename diff_t::Scalar, 1, -1> diagonal;
		int M;
		fdiff_t Kernel;
		template<typename Derived>
		void getKernel(const Eigen::MatrixBase<Derived>&);
		template<typename Derived>
		void KernelModImage(const Eigen::MatrixBase<Derived>&, bool);
	};

	template<typename rank_t, typename diff_t>
	Homology<rank_t, diff_t>::Homology(const Junction<rank_t, diff_t>& J, bool getQ) {
		M = summation(J.rank);
		diff_t In, Out;
		if (J.diffIn.size() == 0 && J.diffOut.size() == 0) {
			In = Eigen::MatrixBase<diff_t>::Zero(M, 1);
			Out = Eigen::MatrixBase<diff_t>::Zero(1, M);
		}
		else if (J.diffOut.size() == 0) {
			Out = Eigen::MatrixBase<diff_t>::Zero(1, M);
			In = J.diffIn;
		}
		else if (J.diffIn.size() == 0) {
			In = Eigen::MatrixBase<diff_t>::Zero(M, 1);
			Out = J.diffOut;
		}
		else {
			Out = J.diffOut;
			In = J.diffIn;
		}

		if constexpr (std::is_integral_v<typename diff_t::Scalar>) 
			getKernel(static_cast<fdiff_t>(Out.template cast<fScalar>()));
		else 
			getKernel(Out);
		if (isZero)
			return;
		fdiff_t In_fScalar = Out_Qi * In.template cast<fScalar>();
		if constexpr (!(std::is_integral_v<typename diff_t::Scalar>)) {
			In = In_fScalar.template cast<typename diff_t::Scalar>();//we cast back if we are using non Z coefficients (F2 etc.). We mustn't do this for Z coefficients due to possible conversion errors.	
			KernelModImage(In, getQ);
		}
		else {
			KernelModImage(In_fScalar, getQ);
		}
	}

	template<typename rank_t, typename diff_t>
	Homology<rank_t, diff_t>::Homology(const Junction<rank_t, diff_t>& J) : Homology(J, 0) {};

	template<typename rank_t, typename diff_t>
	template<typename Derived>
	void Homology<rank_t, diff_t>::getKernel(const Eigen::MatrixBase<Derived>& Out) {

		Smith<Derived, fdiff_t, fdiff_t> OUT(Out, 0, 1);
		Out_Qi = std::move(OUT.Qi);
		Kernel.resize(M,M);
		nonZeroVectors.reserve(M);
		isZero = 1;
		int j = 0;
		for (int i = 0; i < M;i++) {
			if ( (i > OUT.L - 1) || (OUT.diagonal[i] == 0) ) {
				//Remember that Q is invertible so it can't have zero columns
				isZero = 0;
				Kernel.col(j) = OUT.Q.col(i);
				nonZeroVectors.push_back(i);
				j++;
			}
		}
		if (isZero)
			return;
		Kernel.conservativeResize(M, j);
	}

	template<typename rank_t, typename diff_t>
	template<typename Derived>
	void Homology<rank_t, diff_t>::KernelModImage(const Eigen::MatrixBase<Derived>& In, bool getQ) {
		Derived Inner = KeepRow(In, nonZeroVectors);
		auto L = std::min(Inner.rows(), Inner.cols());
	
		Smith<Derived, fdiff_t, fdiff_t> IN(Inner, 1, getQ, 1);
		In_P = std::move(IN.P); 
		if (getQ)
			In_Q = std::move(IN.Q);
		Kernel = Kernel * IN.Pi;
		Generators = Kernel.template cast<typename diff_t::Scalar>();
		auto maxsize = Generators.cols();
		std::vector<typename rank_t::Scalar> groups;
		groups.reserve(maxsize);
		dontModOut.reserve(maxsize);
		isZero = 1;
		diagonal = IN.diagonal.template cast<typename diff_t::Scalar>();
		for (int i = 0; i < maxsize;i++) {
			if (i < L) {
				if (abs(diagonal[i]) == 0) {
					groups.push_back(1);
					isZero = 0;
					dontModOut.push_back(i);
				}
				else if (abs(diagonal[i]) != 1) {
					groups.push_back(static_cast<typename rank_t::Scalar>(abs(diagonal[i])));
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
		if (isZero) {
			Generators.resize(0,0);
			return;
		}
		Groups = Eigen::Map<Eigen::Matrix<typename rank_t::Scalar,1,-1>>(groups.data(), groups.size());
		Generators = KeepCol(Generators, dontModOut);


		//Check for non integer coefficients and replace the Z in Groups with Z/N
		typedef typename diff_t::Scalar coeff;
		if constexpr (is_finite_cyclic<typename diff_t::Scalar>()) {
			constexpr int order = coeff::order;
			for (int i = 0; i < Groups.size(); i++) {
				if (Groups[i] == 1)
					Groups[i] = static_cast<typename rank_t::Scalar>(order);
			}
		}

	}


	template<typename rank_t, typename diff_t>
	template<typename gen_t>
	rank_t Homology<rank_t, diff_t>::basis(const gen_t& generator) const {
		if (isZero)
			return rank_t();
		rank_t basisArray(Groups.size());
		Eigen::Matrix<fScalar,-1,1> element = generator.template cast<fScalar>();
		element = Out_Qi * element;
		element = KeepRow(element, nonZeroVectors);
		element = In_P * element;
		element = KeepRow(element, dontModOut);
		Eigen::Matrix<long, -1, 1> longelement = element.template cast<long>(); //cast long here as it might be large; after that we modulo 

		for (int j = 0; j < Groups.size();j++) {
			typename diff_t::Scalar coeff;
			if (Groups(j) != 1) {
				coeff = static_cast<typename diff_t::Scalar>(longelement(j) % Groups(j));//need to cast incase extra mod is applied by our coefficients. Useless for Z coefficients
				basisArray(j) = (Groups(j) + static_cast<typename rank_t::Scalar>(coeff)) % Groups(j);
				//This is needed due to C++ conventions for % take a symmetric range w.r.t. 0 instead of >=0. So we can't just do element(j)%Groups(j)
			}
			else {
				coeff = static_cast<typename diff_t::Scalar>(longelement(j));//need to cast in case extra mod is applied by our coefficients. Useless for Z coefficients
				basisArray(j) = static_cast<typename rank_t::Scalar>(coeff);
			}
		}
		return basisArray;
	}



	template<typename rank_t, typename diff_t>
	template<typename gen_t>
	gen_t Homology<rank_t, diff_t>::boundary(const gen_t& generator) const {
		Eigen::Matrix<fScalar, -1, 1> element = generator.template cast<fScalar>();
		element = Out_Qi * element;
		element = KeepRow(element, nonZeroVectors);
		element = In_P * element;
		for (const auto& i : dontModOut) {
			if (element[i] != 0) { //the element is not 0 so it has no preimage
				gen_t gen;
				return gen;
			}
		}
		Eigen::Matrix<fScalar, -1, 1> y(In_Q.rows()); //Sy=Px
		y.setZero();
		for (int i = 0; i < y.size(); i++) {
			if (i < diagonal.size() && diagonal[i] != 0)
				y[i] = static_cast<fScalar>(element[i] / diagonal[i]);
		}
		element = In_Q * y;
		return element.template cast<typename gen_t::Scalar>();
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
	void mod_normalize(Eigen::Matrix<T, 1, -1>& element, const Eigen::Matrix<T, 1, -1>& group) {
		for (int i = 0; i < group.size(); i++) {
			if (group[i] != 1)
				element[i] = (element[i] + group[i]) % group[i];
		}
	}


	///Normalizes the element in a group to have minimal signs amongst its multiples
	template<typename T>
	void normalize(Eigen::Matrix<T, 1, -1>& element, const Eigen::Matrix<T, 1, -1>& group) {
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
			typename T::Scalar q = 0;
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
		for (typename T::Scalar i = 1; i < ob; i++) {
			if (equals(a, i * b, group))
				return i;
		}
		return 0;
	}

}
