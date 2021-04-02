#pragma once
#include "../Homology.hpp"

///@file
///@brief Contains the Homology class and algorithms.

namespace mackey {

	template<typename rank_t, typename diff_t>
	Homology<rank_t, diff_t>::Homology(const Junction<rank_t, diff_t>& J, bool getQ) {
		M = summation(J.rank);
		diff_t_C In, Out;
		if (J.diffIn.size() == 0) {
			In.resize(M, 1);
			In.setZero();
		}
		else
			In= J.diffIn.template cast<HScalar>();
		if (J.diffOut.size() == 0) {
			Out.resize(1,M);
			Out.setZero();
		}
		else 
			Out = J.diffOut.template cast<HScalar>();
		auto Kernel=getKernel(Out);
		if (isZero)
			return;
		KernelModImage(In, Kernel, getQ);
	}

	template<typename rank_t, typename diff_t>
	typename Homology<rank_t, diff_t>::diff_t_C Homology<rank_t, diff_t>::getKernel(diff_t_C& Out) {
		Smith_Normal_Form<diff_t_C, diff_t_R, diff_t_C> OUT(Out, 0, 1, 0);
		diff_t_C Kernel(M,M);
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
	void Homology<rank_t, diff_t>::KernelModImage(diff_t_C& In, diff_t_C& Kernel, bool getQ) {
		if constexpr (SFINAE::is_Sparse<diff_t>::value)
			In = (Out_Qi * In).pruned();
		else
			In = Out_Qi * In;
		auto L = std::min(In.rows(), In.cols());
		Smith_Normal_Form<diff_t_C, diff_t_R, diff_t_C> IN(In, 1, getQ, 1);

		if constexpr (SFINAE::is_Sparse<diff_t>::value)
			Generators= (Kernel * IN.Pi).pruned();
		else
			Generators = Kernel * IN.Pi;
		auto maxsize = Generators.cols();
		std::vector<scalar_t<rank_t>> groups;
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
					groups.push_back(static_cast<scalar_t<rank_t>>(abs(diagonal[i])));
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
		rank_t G=Eigen::Map<rank_t>(groups.data(), groups.size());
		Groups = AbelianGroup(G);
		Generators = KeepCol(Generators, dontModOut);


		//Check for non integer coefficients and replace the Z in Groups with Z/N
		if constexpr (SFINAE::is_finite_cyclic<scalar_t<diff_t>>::value) {
			constexpr int order = scalar_t<diff_t>::order;
			for (int i = 0; i < Groups.number_of_summands(); i++) {
				if (Groups[i] == 1)
					Groups[i] = static_cast<scalar_t<rank_t>>(order);
			}
		}
	}

	template<typename rank_t, typename diff_t>
	rank_t Homology<rank_t, diff_t>::basis(const gen_t& generator) const {
		if (isZero)
			return rank_t();
		rank_t basisArray(Groups.number_of_summands());
		gen_t element = Out_Qi * generator;
		element = In_P_reduced * element;
		for (int j = 0; j < Groups.number_of_summands();j++) {
			if (Groups[j] != 1)
				basisArray[j] = static_cast<scalar_t<rank_t>>((Groups[j] + (int64_t)element[j] % Groups[j]) % Groups[j]);
			else
				basisArray[j] = static_cast<scalar_t<rank_t>>(element[j]);
		}
		return basisArray;
	}

	template<typename rank_t, typename diff_t>
	typename Homology<rank_t, diff_t>::gen_t Homology<rank_t, diff_t>::boundary(const gen_t& generator) const {
		gen_t element = In_P_full * Out_Qi * generator;
		for (const auto& i : dontModOut) {
			if ((diagonal[i] != 0 && (int64_t)element[i]%(int64_t)diagonal[i] != 0) || (diagonal[i]==0 && element[i]!=0))  //the element is not 0 so it has no preimage
				return gen_t();
		}
		gen_t y(In_Q.rows()); //Sy=Px
		y.setZero();
		for (typename gen_t::StorageIndex i = 0; i < y.size(); i++) {
			if (i < diagonal.size() && diagonal[i] != 0)
				y[i] = static_cast<HScalar>(element[i] / diagonal[i]);
		}
		return In_Q * y;
	}
}
