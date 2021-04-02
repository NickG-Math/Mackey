#pragma once
#include "../Green.hpp"

///@file
///@brief Contains the class and methods for multiplying generators.

namespace mackey {

	template<typename group_t>
	Green<group_t>::Green() { isZero = 1; }

	template<typename group_t>
	Green<group_t> ::Green(const chains_t<group_t>& C, const chains_t<group_t>& D, int level, int degreeC, int degreeD)
	{
		internal::GreenCompute<group_t> Comp(C, D, level, degreeC, degreeD);
		*this = std::move(Comp.G);
	}

	template<typename group_t>
	template<typename T, typename S>
	auto Green<group_t>::getNormalBasis(int i, const S& b) const {
		if (isZero)
			return rank_t();
		return getNormalBasis(basisElement<rank_t>(first_number_selections, i), b);
	}

	template<typename group_t>
	template<typename T, typename S>
	auto Green<group_t>::getNormalBasis(const T& a, int j) const {
		if (isZero)
			return rank_t();
		return getNormalBasis(a, basisElement<rank_t>(second_number_selections, j));
	}

	template<typename group_t>
	template<typename T, typename S>
	auto Green<group_t>::getNormalBasis(int i, int j) const {
		if (isZero)
			return rank_t();
		auto b = basis[select(i, j)];
		Groups.normalize(b);
		return b;
	}

	template<typename group_t>
	int Green<group_t>::select(int select_first, int select_second) const {
		if (select_first >= first_number_selections)
			select_first = 0;
		if (select_second >= second_number_selections)
			select_second = 0;
		return select_second + select_first * second_number_selections;
	}

	template<typename group_t>
	template<typename T, typename S>
	auto Green<group_t>::getNormalBasis(const T& c, const S& d) const {
		if (isZero)
			return rank_t();
		rank_t result(basis[0].size());
		result.setZero();
		for (int i = 0; i < c.size(); i++) {
			if (c[i] == 0)
				continue;
			for (int j = 0; j < d.size(); j++) {
				if (d[j] == 0)
					continue;
				result += c[i] * d[j] * basis[select(i, j)];
			}
		}
		Groups.normalize(result);
		return result;
	}

	template<typename group_t>
	bool Green<group_t>::operator==(const Green<group_t>& other) const {
		return isZero == other.isZero && Groups.size() == other.Groups.size() && Groups == other.Groups && basis.size() == other.basis.size() && basis == other.basis && boxID == other.boxID;
	}


	template<typename group_t>
	typename group_t::rank_t multiply(const chains_t<group_t>& C, const chains_t<group_t>& D, int level, int degreeC, int degreeD, int selectC, int selectD){
		Green<group_t> G(C, D, level, degreeC, degreeD);
		return G.getNormalBasis(selectC, selectD);
	}

	namespace internal {

		///Computes the homology at given level, the generators and their restrictions
		template<typename group_t>
		ChainsLevelGen<group_t>::ChainsLevelGen(const Chains<rank_t, diff_t>& C, int level, int degree) {
			Junction<rank_t, diff_t> J(C, degree);
			Junction<rank_t, diff_t> J_Level = transfer<group_t>(J, level);
			Homology<rank_t, diff_t> Homol(J_Level);
			isZero = Homol.isZero;
			if (isZero)
				return;
			Gens = std::move(Homol.Generators);
			rank_level = std::move(J_Level.rank);
		}

		template<typename group_t>
		ChainsLevelGen<group_t>::ChainsLevelGen(const Chains<rank_t, diff_t>& C, int level, int degree, int selection) : ChainsLevelGen(C, level, degree) {
			if (Gens.cols() > 1)
				gen = Gens.col(selection);
			else
				gen = Gens.col(0);
			res_gen = restriction(gen, rank_level, C.rank[degree]);
		}

		///Computes the product of generators at the bottom level
		template<typename rank_t, typename diff_t>
		ProductGen<rank_t, diff_t>::ProductGen(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int degreeC, int degreeD) : C(C), D(D), degreeC(degreeC), degreeD(degreeD) {}

		template<typename rank_t, typename diff_t>
		void ProductGen<rank_t, diff_t>::pad(const std::vector<int64_t>& detailedrank) {
			padleft = padright = 0;
			for (int i = 0; i <= degreeC - 1; i++)
				padleft += detailedrank[i];
			for (int i = degreeC + 1; i <= degreeC + degreeD; i++)
				padright += detailedrank[i];
		}
		template<typename rank_t, typename diff_t>
		void ProductGen<rank_t, diff_t>::multiply(gen_t<rank_t, diff_t> gen1, gen_t<rank_t, diff_t> gen2) {
			gen_t<rank_t, diff_t> leftConvProduct(gen1.size() * gen2.size());
			for (int i = 0; i < gen2.size(); i++)
				leftConvProduct.segment(i * gen1.size(), gen1.size()) = gen1 * gen2(i);
			ChangeBasis<int> permute(C.rank[degreeC], D.rank[degreeD], 1, 0);
			const Eigen::PermutationMatrix<-1, -1, int> left_to_canon(Eigen::Map<const Eigen::Matrix<int, -1, 1>>(permute.left_to_canon.data(), permute.left_to_canon.size()));
			gen_t<rank_t, diff_t> canonProduct = left_to_canon * leftConvProduct;
			product.resize(padleft + canonProduct.size() + padright);
			product.setZero();
			product.segment(padleft, canonProduct.size()) = canonProduct;
		}


		template<typename group_t>
		GreenCompute<group_t>::GreenCompute(const chains_t<group_t>& C, const chains_t<group_t>& D, int level, int degreeC, int degreeD)
			: C(C), D(D), level(level), degreeC(degreeC), degreeD(degreeD), C_Variables(ChainsLevelGen<group_t>(C, level, degreeC)), D_Variables(ChainsLevelGen<group_t>(D, level, degreeD)), Restricted(ProductGen(C, D, degreeC, degreeD))
		{
			getGens();
			if (G.isZero)
				return;
			box();
			if (G.isZero)
				return;
			compute();
		}

		template<typename group_t>
		void GreenCompute<group_t>::getGens() {
			GenC = std::move(C_Variables.Gens);
			if (C_Variables.isZero)
				return;
			GenD = std::move(D_Variables.Gens);
			G.first_number_selections = GenC.cols();
			G.second_number_selections = GenD.cols();
			G.isZero = C_Variables.isZero & D_Variables.isZero;
		}

		template<typename group_t>
		void GreenCompute<group_t>::box() {
			Tensor<Junction<rank_t, diff_t>, std::vector<int64_t>> Boxed(C, D, degreeC + degreeD);
			rank_bottom = Boxed.tensor().rank;
			Restricted.pad(Boxed.optional());
			IDGeneratorCompute<group_t> ID(level, Boxed.tensor());
			H_Level = std::move(ID.H_level);
			rank_level = std::move(ID.rank_level);
			G.boxID = std::move(ID.ID);
			G.isZero = H_Level.isZero;
		}

		template<typename group_t>
		void GreenCompute<group_t>::compute() {
			gen_t<rank_t, diff_t> genC, genD;
			auto combinations = GenC.cols() * GenD.cols();
			G.basis.reserve(combinations);
			for (int i = 0; i < GenC.cols(); i++) {
				for (int j = 0; j < GenD.cols(); j++) {
					genC = GenC.col(i);
					genD = GenD.col(j);
					auto resGenC = restriction(genC, C_Variables.rank_level, C.rank[degreeC]);
					auto resGenD = restriction(genD, D_Variables.rank_level, D.rank[degreeD]);
					Restricted.multiply(resGenC, resGenD);
					product = invRes(Restricted.product, rank_bottom, rank_level);
					rank_t basis = H_Level.basis(product);
					G.basis.push_back(basis);
				}
			}
			G.Groups = std::move(H_Level.Groups);
		}
	}


}
