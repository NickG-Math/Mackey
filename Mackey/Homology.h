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
<<<<<<< HEAD
=======
		
//////////////////////////////////////////////////
///The floating point type for Homology computations

///Should be float to maximize performance, or potentially double for groups of order not a prime power.
/////////////////////////////////////////////////
		typedef float fScalar;
		typedef Eigen::Matrix<fScalar, -1, -1> fdiff_t;	///<The type of our casted differentials


>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
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
<<<<<<< HEAD
			Out = Eigen::MatrixBase<diff_t>::Zero(1, M);
			In = J.diffIn;
=======
			Out = Eigen::MatrixBase<fdiff_t>::Zero(1, M);
			In = J.diffIn.template cast<fScalar>();
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
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
		if constexpr (! (std::is_integral_v<typename diff_t::Scalar>)) {
			In = In_fScalar.template cast<typename diff_t::Scalar>();//we cast back if we are using non Z coefficients (F2 etc.). We mustn't do this for Z coefficients due to possible conversion errors.	
			KernelModImage(In, getQ);
		}
		else 
			KernelModImage(In_fScalar, getQ);
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
<<<<<<< HEAD
			if ( (i > OUT.L - 1) || (OUT.diagonal[i] == 0) ) {
=======
			if (((i > OUT.L - 1) || (static_cast<typename diff_t::Scalar>(OUT.S(i, i)) == 0))) {
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
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
		Smith<Derived, fdiff_t, fdiff_t> IN(Inner, 1, getQ);
		In_P = std::move(IN.P); 
		if (getQ)
			In_Q = std::move(IN.Q);
		Kernel = Kernel * IN.Pi;
		Generators = Kernel.template cast<typename diff_t::Scalar>();
		auto maxsize = Generators.cols();
<<<<<<< HEAD
		std::vector<typename rank_t::Scalar> groups;
		groups.reserve(maxsize);
=======
		Groups.resize(maxsize);
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
		dontModOut.reserve(maxsize);
		isZero = 1;
		diagonal = IN.diagonal.template cast<typename diff_t::Scalar>();
		for (int i = 0; i < maxsize;i++) {
<<<<<<< HEAD
			if (i < L) {
				if (abs(diagonal[i]) == 0) {
					groups.push_back(1);
					isZero = 0;
					dontModOut.push_back(i);
				}
				else if (abs(diagonal[i]) != 1) {
					groups.push_back(static_cast<typename rank_t::Scalar>(abs(diagonal[i])));
=======
			if (i < In.rows() && i < In.cols()) {
				auto diagonalentry = static_cast<typename diff_t::Scalar>(abs(IN.S(i, i)));
				if (diagonalentry == 0) {
					Groups(i) = 1;
					isZero = 0;
					dontModOut.push_back(i);
				}
				else if (diagonalentry != 1) {
					Groups(i) = static_cast<typename rank_t::Scalar>(diagonalentry);
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
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
<<<<<<< HEAD
		Groups = Eigen::Map<Eigen::Matrix<typename rank_t::Scalar,1,-1>>(groups.data(), groups.size());
=======
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
		Generators = KeepCol(Generators, dontModOut);
	}


	template<typename rank_t, typename diff_t>
	template<typename gen_t>
	rank_t Homology<rank_t, diff_t>::basis(const gen_t& generator) const {
		if (isZero) {
			rank_t basisArray;
			return basisArray;
		}
		rank_t basisArray(Groups.size());
		Eigen::Matrix<fScalar,-1,1> element = generator.template cast<fScalar>();
		element = Out_Qi * element;
		element = KeepRow(element, nonZeroVectors);
		element = In_P * element;
		element = KeepRow(element, dontModOut);
<<<<<<< HEAD
		Eigen::Matrix<long, -1, 1> longelement = element.template cast<long>(); //cast long here as it might be large; after that we modulo 
=======
		Eigen::Matrix<typename diff_t::Scalar, -1, 1> castelement = element.template cast<typename diff_t::Scalar>();
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
		for (int j = 0; j < Groups.size();j++) {
			typename diff_t::Scalar coeff;
			if (Groups(j) != 1) {
<<<<<<< HEAD
				coeff = static_cast<typename diff_t::Scalar>(longelement(j) % Groups(j));//need to cast incase extra mod is applied by our coefficients. Useless for Z coefficients
				basisArray(j) = (Groups(j) + static_cast<typename rank_t::Scalar>(coeff)) % Groups(j);
=======
				basisArray(j) = (Groups(j) + static_cast<typename rank_t::Scalar>(castelement(j)) % Groups(j)) % Groups(j);
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
				//This is needed due to C++ conventions for % take a symmetric range w.r.t. 0 instead of >=0. So we can't just do element(j)%Groups(j)
				//Also we can't use our own mod since this is mod of floats over int (although we could just define mod to just cast to int.
			}
			else {
<<<<<<< HEAD
				coeff = static_cast<typename diff_t::Scalar>(longelement(j));//need to cast incase extra mod is applied by our coefficients. Useless for Z coefficients
				basisArray(j) = static_cast<typename rank_t::Scalar>(coeff);
=======
				basisArray(j) = static_cast<typename rank_t::Scalar>(castelement(j));
>>>>>>> df2e02523d9b0cf2be901620a7daf93ba97830a8
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



///Normalizes the element in given abelian group so that element[i]=1 if it generates group[i]
	template<typename T, typename S>
	T normalize(const T& element, const S& group) {
		T normalized = element;
		for (int i = 0; i < element.size(); i++) {
			if (group[i] != 1 && element[i] != 0) {
				normalized[i] = std::gcd(element[i], group[i]);
			}
			else {
				normalized[i] = abs(element[i]);
			}
		}
		return normalized;
	}


	////////////////////////////////////////////////////
///Finds the order of element in given abelian group.

///T,S are vectors or 1d Eigen matrices
///It's assumed that element has been normalized i.e. element[i]=1 if it generates group[i]
///////////////////////////////////////////////////
	template<typename T, typename S>
	int order(const T& element, const S& group) {

		std::vector<int> fractions;
		fractions.reserve(group.size());

		for (int i = 0; i < group.size(); i++) {
			if (group[i] == 1 && element[i] != 0) {
				return 1; //infinite order
			}
			else if (element[i] != 0) {
				fractions.push_back(group[i] / element[i]);
			}
		}
		return lcm(fractions);
	}

	///If the given element is 0,...,0,?,0,...,0 returns the index of ?. Otherwise returns -1
	template<typename T>
	int isBasisElement(const T& element) { // >=0 if it's 0,...,0,?,0,...,0
		int nonZero = 0;
		int pos = 0;
		for (int i = 0; i < element.size(); i++) {
			if (element[i] != 0) {
				nonZero++;
				pos = i;
				if (nonZero>1)
					return -1; 
			}
		}
		return pos;
	}

	///Given an array of elements, returns the index of the unique element of the form 0,...,0,?,0,...,0. If not unique then returns -1
	template<typename T>
	int findBasisElement(const T& elements) { // >=0 if it's 0,...,0,?,0,...,0
		int counter = 0;
		int pos = 0;
		for (std::vector<int>::size_type i=0; i<elements.size(); i++){
			if (isBasisElement(elements[i]) != -1) {
				counter++;
				pos = i;
				if (counter > 1)
					return -1;
			}
		}
		if (counter==0)
			return -1;
		return pos;
	}
}
