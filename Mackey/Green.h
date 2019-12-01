#pragma once
#include "Levels.h"
#include "Box.h"
#include "Identify.h"
#include <numeric> //for gcd

///@file
///@brief Contains the class and methods for multiplying generators.

namespace Mackey {

	/// The result of multiplying generators in a Green functor
	template<typename rank_t, typename diff_t>
	class Green {
	public:
		rank_t Groups;///<The homology group the product lives in.
		bool isZero;///<1 if the the homology group the product lives in is 0.

		//////////////////////////////////////////////////////////////
		///Expresses the product of two generators as a linear combination of the generators of the box product

		std::vector<rank_t> basis;		///<For each selection of generators we have an array expressing the product as a linear combination of generators in the degree of the product
		std::vector<rank_t> normalBasis; ///<Same as basis, now normalized.

		int first_number_selections;///<The number of selections of generators for the first factor
		int second_number_selections;///<The number of selections of generators for the second factor

		IDGenerators<rank_t, diff_t> boxID;///<Identifies the generators of the homology group the product lives in.
		
		///Default constructor
		Green() { isZero = 1; };

		///Computes the product of generators given chains, degrees and level
		Green(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int, int, int);

		///Selects the desired generator
		int select(int, int) const;

	};


	template<typename rank_t, typename diff_t>
	int Green<rank_t, diff_t>::select(int select_first, int select_second) const {
		if (select_first >= first_number_selections)
			select_first = 0;
		if (select_second >= second_number_selections)
			select_second = 0;
		return select_second + select_first * second_number_selections;
	}


	namespace internal {
		
		///Computes the homology at given level, the generators and their restrictions
		template<typename rank_t, typename diff_t>
		class ChainsLevelGen {
			typedef Eigen::Matrix<typename diff_t::Scalar, -1, 1> gen_t;
			diff_t Gens;
			gen_t gen, res_gen;
			rank_t rank_level;
			bool isZero;
			ChainsLevelGen(const Chains<rank_t, diff_t>& C, int level, int degree) {
				Junction<rank_t, diff_t> J(C, degree);
				Junction<rank_t, diff_t> J_Level = transfer(J, level);
				Homology<rank_t, diff_t> Homol(J_Level);
				isZero = Homol.isZero;
				if (isZero)
					return;
				Gens = std::move(Homol.Generators);
				rank_level = std::move(J_Level.rank);
			}
			ChainsLevelGen(const Chains<rank_t, diff_t>& C, int level, int degree, int selection)
				: ChainsLevelGen(C, level, degree)
			{
				if (Gens.cols() > 1)
					gen = Gens.col(selection);
				else
					gen = Gens.col(0);
				res_gen = restriction(gen, rank_level, C.rank[degree]);
			}

			template<typename s_rank_t, typename s_diff_t>
			friend class GreenCompute;
			
			template<typename s_rank_t, typename s_diff_t>
			friend class MasseyCompute;

		};


		template<typename rank_t, typename diff_t>
		class ProductGen {
		public:
			typedef Eigen::Matrix<typename diff_t::Scalar, -1, 1> gen_t;
			gen_t product;
			int padright, padleft;
			const Chains<rank_t, diff_t>& C, D;
			int degreeC, degreeD;
			ProductGen(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int degreeC, int degreeD) : C(C), D(D), degreeC(degreeC), degreeD(degreeD) {}

			void pad(const JunctionBox<rank_t, diff_t>& Boxed) {
				pad(Boxed.detailedrank);
			}

			void pad(const std::vector<rank_t>& detailedrank) {
				padleft = padright = 0;
				for (int i = 0; i <= degreeC - 1; i++) {
					padleft += summation(detailedrank[i]);
				}
				for (int i = degreeC + 1; i <= degreeC + degreeD; i++) {
					padright += summation(detailedrank[i]);
				}
			}


			void multiply(gen_t gen1, gen_t gen2) {
				gen_t leftConvProduct(gen1.size() * gen2.size());
				for (int i = 0; i < gen2.size(); i++) {
					leftConvProduct.segment(i * gen1.size(), gen1.size()) = gen1 * gen2(i);
				}
				ChangeBasis<rank_t> permute(C.rank[degreeC], D.rank[degreeD]);
				gen_t canonProduct = permute.LefttoCanon.inverse() * leftConvProduct;
				product.resize(padleft + canonProduct.size() + padright);
				product.setZero();
				product.segment(padleft, canonProduct.size()) = canonProduct;
			}

		};

		/// Computes the product of generators in a Green functor
		template<typename rank_t, typename diff_t>
		class GreenCompute {

			Green<rank_t, diff_t> G; ///<The answer of the computation

			typedef Eigen::Matrix<typename diff_t::Scalar, -1, 1> gen_t;
			gen_t product;


			const Chains<rank_t, diff_t>& C, D;
			diff_t GenC, GenD;
			rank_t rank_bottom, rank_level;
			const int level, degreeC, degreeD;
			Homology<rank_t, diff_t> H_Level;
			ChainsLevelGen<rank_t, diff_t> C_Variables, D_Variables;
			ProductGen<rank_t, diff_t> Restricted;

			///Computes the product of generators
			GreenCompute(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int level, int degreeC, int degreeD);


			/// Computes the generators from their degrees and the selections
			void getGens();
			/// Computes the tensor product of Chains C,D
			void box();
			/// Computes the product of generators genC, genD by combining the results of getGens and box
			void compute();

			template<typename s_rank, typename s_diff>
			friend class Mackey::Green;
		};


		template<typename rank_t, typename diff_t>
		GreenCompute<rank_t, diff_t> ::GreenCompute(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int level, int degreeC, int degreeD)
			: C(C), D(D), level(level), degreeC(degreeC), degreeD(degreeD), C_Variables(ChainsLevelGen(C, level, degreeC)), D_Variables(ChainsLevelGen(D, level, degreeD)), Restricted(ProductGen(C, D, degreeC, degreeD))
		{
			getGens();
			if (G.isZero)
				return;
			box();
			if (G.isZero)
				return;
			compute();
		}

		template<typename rank_t, typename diff_t>
		void GreenCompute<rank_t, diff_t>::getGens() {
			GenC = std::move(C_Variables.Gens);
			if (C_Variables.isZero)
				return;
			GenD = std::move(D_Variables.Gens);
			G.first_number_selections = GenC.cols();
			G.second_number_selections = GenD.cols();
			G.isZero = C_Variables.isZero & D_Variables.isZero;
		}

		template<typename rank_t, typename diff_t>
		void GreenCompute<rank_t, diff_t>::box() {
			JunctionBox<rank_t, diff_t> Boxed(C, D, degreeC + degreeD);
			rank_bottom = Boxed.rank;
			Restricted.pad(Boxed);
			IDGeneratorCompute<rank_t, diff_t> ID(level, Boxed);
			H_Level = std::move(ID.H_level);
			rank_level = std::move(ID.rank_level);
			G.boxID = std::move(ID.ID);
			G.isZero = H_Level.isZero;
		}

		template<typename rank_t, typename diff_t>
		void GreenCompute<rank_t, diff_t>::compute() {

			typedef Eigen::Matrix<typename diff_t::Scalar, -1, 1> gen_t;
			gen_t genC, genD;
			auto combinations = GenC.cols() * GenD.cols();
			G.basis.reserve(combinations);
			G.normalBasis.reserve(combinations);
			for (int i = 0; i < GenC.cols(); i++) {
				for (int j = 0; j < GenD.cols(); j++) {
					genC = GenC.col(i);
					genD = GenD.col(j);
					gen_t resGenC = restriction(genC, C_Variables.rank_level, C.rank[degreeC]);
					gen_t resGenD = restriction(genD, D_Variables.rank_level, D.rank[degreeD]);
					Restricted.multiply(resGenC, resGenD);
					product = invRes(Restricted.product, rank_bottom, rank_level);
					rank_t basis = H_Level.basis(product);

					G.basis.push_back(basis);

					auto normalBasis = normalize(basis, H_Level.Groups);
					G.normalBasis.push_back(normalBasis);
				}
			}
			G.Groups = std::move(H_Level.Groups);
		}

	}

	template<typename rank_t, typename diff_t>
	Green<rank_t, diff_t> ::Green(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int level, int degreeC, int degreeD)
	{
		internal::GreenCompute Comp(C, D, level, degreeC, degreeD);
		*this = std::move(Comp.G);
	}


	template<typename rank_t, typename diff_t>
	rank_t multiply(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int level, int degreeC, int degreeD, int selectC, int selectD)
	{
		Green G(C, D, level, degreeC, degreeD);
		return G.normalBasis[G.select(selectC, selectD)];
	}

	template<typename rank_t, typename diff_t>
	rank_t multiply(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, int level, int degreeC, int degreeD)
	{
		Green G(C, D, level, degreeC, degreeD);
		return G.normalBasis[0];
	}


}
