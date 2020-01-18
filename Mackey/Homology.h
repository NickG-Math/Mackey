#pragma once
#include "Chains.h"
#include "General.h"
#include "Smith.h"
#include "Z_n.h"
#include <numeric>

///@file
///@brief Contains the Homology class and algorithms.

namespace Mackey {
	
	namespace {
		//to avoid too many std conditionals

		template<typename, typename = void>
		struct Quicktypedefs;

		template<typename T>
		struct Quicktypedefs<T, typename std::enable_if_t<SFINAE::is_Dense<T>::value && SFINAE::is_finite_cyclic<Scalar_t<T>>::value>> {
			typedef Scalar_t<T> Scalar;
			typedef Eigen::Matrix<Scalar, -1, -1> col;
			typedef Eigen::Matrix<Scalar, -1, -1, 1> row;
		};

		template<typename T>
		struct Quicktypedefs<T, typename std::enable_if_t<SFINAE::is_Sparse<T>::value && SFINAE::is_finite_cyclic<Scalar_t<T>>::value>> {
			typedef Scalar_t<T> Scalar;
			typedef Eigen::SparseMatrix<Scalar, 0> col;
			typedef Eigen::SparseMatrix<Scalar, 1> row;
		};

		template<typename T>
		struct Quicktypedefs<T, typename std::enable_if_t<SFINAE::is_Dense<T>::value && !SFINAE::is_finite_cyclic<Scalar_t<T>>::value>> {
			typedef float Scalar;
			typedef Eigen::Matrix<Scalar, -1, -1> col;
			typedef Eigen::Matrix<Scalar, -1, -1, 1> row;
		};

		template<typename T>
		struct Quicktypedefs<T, typename std::enable_if_t<SFINAE::is_Sparse<T>::value && !SFINAE::is_finite_cyclic<Scalar_t<T>>::value>> {
			typedef long Scalar;
			typedef Eigen::SparseMatrix<Scalar, 0, typename T::StorageIndex> col;
			typedef Eigen::SparseMatrix<Scalar, 1, typename T::StorageIndex> row;
		};


	}

	///The Homology of a Junction
	template<typename rank_t, typename diff_t>
	class Homology {

		///<Type used for Smith + product of matrices, either for performance (float/double) or to avoid overflow (long).
		typedef typename Quicktypedefs<diff_t>::Scalar fScalar;
		typedef typename Quicktypedefs<diff_t>::col fdiff_t_C;	///<Column Major
		typedef typename Quicktypedefs<diff_t>::row fdiff_t_R;	///<Row Major
	public:

		///The type of our matrix of generators
		typedef fdiff_t_C Gens_t;
		///The dense type of our generators (a column in the generator matrix, always dense for convenience)
		typedef col_t<fdiff_t_C> gen_t;


		rank_t Groups;///<Encodes the homology groups as follows: Groups=[1,2,3] means homology Z+Z/2+Z/3
		Gens_t Generators;///<Encodes the generators homology groups as follows: The i-th column corresponds to the generator for Groups[i]
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
		std::vector<typename diff_t::StorageIndex> dontModOut;
		fdiff_t_C In_Q;
		fdiff_t_R Out_Qi, In_P_full, In_P_reduced;
		row_t<fdiff_t_C> diagonal;
		typename diff_t::StorageIndex M;
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
		std::vector<typename diff_t::StorageIndex> nonZeroVectors;
		nonZeroVectors.reserve(M);
		isZero = 1;
		typename diff_t::StorageIndex j = 0;
		for (typename diff_t::StorageIndex i = 0; i < M;i++) {
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
		if constexpr (SFINAE::is_Sparse<diff_t>::value)
			In = (Out_Qi * In).pruned().eval();
		else
			In = Out_Qi * In;
		auto L = std::min(In.rows(), In.cols());
		auto IN=diagonalize<fdiff_t_C, fdiff_t_R, fdiff_t_C>(In, 1, getQ, 1);
		if constexpr (SFINAE::is_Sparse<diff_t>::value)
			Generators= (Kernel * IN.Pi).pruned();
		else
			Generators = Kernel * IN.Pi;
		auto maxsize = Generators.cols();
		std::vector<Scalar_t<rank_t>> groups;
		groups.reserve(maxsize);
		dontModOut.reserve(maxsize);
		isZero = 1;
		diagonal = std::move(IN.diagonal);
		for (typename diff_t::StorageIndex i = 0; i < maxsize;i++) {
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
		Generators = KeepCol(Generators, dontModOut);
		Groups = Eigen::Map<rank_t>(groups.data(), groups.size());


		//Check for non integer coefficients and replace the Z in Groups with Z/N
		if constexpr (SFINAE::is_finite_cyclic<Scalar_t<diff_t>>::value) {
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
		gen_t element = Out_Qi * generator;
		element = (In_P_reduced * element).eval();
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
			if ((diagonal[i] != 0 && (long)element[i]%(long)diagonal[i] != 0) || (diagonal[i]==0 && element[i]!=0))  //the element is not 0 so it has no preimage
				return gen_t();
		}
		gen_t y(In_Q.rows()); //Sy=Px
		y.setZero();
		for (typename gen_t::StorageIndex i = 0; i < y.size(); i++) {
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
///Expresses a as a linear combination of elements b[i] living in a finitely generated abelian group. If a is not in their span, returns empty matrix.
	template<typename T>
	auto span(const T& a, const std::vector<T>& b, const T& group) {
		int count = 0;
		for (int i = 0; i < group.size(); i++) {
			if (group[i] != 1)
				count++;
		}
		mat_t<T> relation(group.size(), count + b.size());
		relation.setZero();
		int j = 0;
		for (int i = 0; i < group.size(); i++) {
			if (group[i] != 1) {
				relation(j, j) = group[i];
				j++;
			}
		}
		for (int i = 0; i < b.size(); i++)
			relation.col(i + count) = b[i];
		Junction<T, mat_t<T>> J;
		J.rank = T::Constant(1, (int)group.size());
		J.diffIn = relation;
		Homology<T, mat_t<T>> H(J, 1);
		gen_t<T, mat_t<T>> cast_a = a.template cast<Scalar_t<gen_t<T, mat_t<T>>>>();
		auto v = H.boundary(cast_a);
		if (v.size()==0)
			return v;
		decltype(v) z(b.size());
		for (int i = 0; i < b.size(); i++)
			z[i] = v[i+count];
		return z;
	}

	////////////////////////////////////////////////////
///Finds if an element is in the span of other given elements in a finitely generated abelian group
	template<typename T>
	bool inSpan(const T& a, const std::vector<T>& b, const T& group) {
		if (span(a, b, group).size() != 0)
			return 1;
		return 0;
	}


	////////////////////////////////////////////////////
///Returns k if a=k*b in group. Returns 0 if a is not a multiple of b (so make sure a!=0)
	template<typename T>
	int isMultiple(const T& a, const T& b, const T& group) {
		auto u = span(a, { b }, group);
		if (u.size() == 0)
			return 0;
		return (int)u[0];
	}
}
