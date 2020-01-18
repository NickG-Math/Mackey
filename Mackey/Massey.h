#pragma once
#include "Green.h"

///@file
///@brief Contains the class and methods for Massey Products.

namespace Mackey {

	///Stores Massey products and their indeterminacy
	template<typename rank_t, typename diff_t>
	class Massey {
	public:
		rank_t basis; ///<Expresses the Massey product as a linear combination of the generators of the box product.
		rank_t normalBasis; ///<Same as basis, but normalized
		rank_t Groups; ///<The homology group the product lives in.
		std::array<rank_t, 2> indeterminacy; ///< The indeterminacy of the Massey product, stored as two bases
		IDGenerators<rank_t> boxID; ///< Identification data for the generators of the Massey product
		bool exists; ///< If the Massey product exists
		bool noIndeterminacy; ///<If the Massey product has no indeterminacy
		bool isZero; ///<If the Massey product is 0.
	
		///Default constructor
		Massey() { exists = 0; isZero = 1; }

		///Compute Massey product given level, generator degrees and selections
		Massey(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int, int, int, int, int, int, int);

		///Compute Massey product given level, generator degrees and default 0,0,0 selections
		Massey(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int, int, int, int);

	};

	namespace internal 
	{

		///Forms the permutation allowing us to go from C box (D box E) to (C box D) box E.
		template<typename rank_t>
		Eigen::PermutationMatrix<-1, -1, int> box_permutation(const std::vector<rank_t>& C, const std::vector<rank_t>& D, const std::vector<rank_t>& E, int degree) {
			auto Cmax = std::min((int)C.size() - 1, degree);
			auto Dmax = std::min((int)D.size() - 1, degree);
			auto Emax = std::min((int)E.size() - 1, degree);
			std::vector<int> CD_E;
			for (int k = Emax; k >= 0; k--) {
				for (int j = Dmax; j >= 0; j--) {
					for (int i = 0; i <= Cmax; i++) {
						if (i + j + k == degree) {
							int rank = summation(C[i]) * summation(D[j]) * summation(E[k]);
							for (int u = 0; u < rank; u++) {
								CD_E.push_back(i + j * (Cmax + 1) + k * (Cmax + 1) * (Dmax + 1) + u * (Cmax + 1) * (Dmax + 1) * (Emax + 1)); //k,i,j, pushing back as many times as the rank
							}
						}
					}
				}
			}
			std::vector<int> C_DE;
			for (int i = 0; i <= Cmax; i++) {
				for (int j = 0; j <= Dmax; j++) {
					for (int k = 0; k <= Emax; k++) {
						if (i + j + k == degree) {
							int rank = summation(C[i]) * summation(D[j]) * summation(E[k]);
							for (int u = 0; u < rank; u++) {
								C_DE.push_back(i + j * (Cmax + 1) + k * (Cmax + 1) * (Dmax + 1) + u * (Cmax + 1) * (Dmax + 1) * (Emax + 1)); //i,j,k
							}
						}
					}
				}
			}
			auto perm = changebasis<int>(CD_E,C_DE);
			return Eigen::PermutationMatrix<-1, -1, int>(Eigen::Map<Eigen::Matrix<int, -1, 1>>(perm.data(), perm.size()));
		}

		///Computes Massey products and their indeterminacy
		template<typename rank_t, typename diff_t>
		class MasseyCompute {
			Massey<rank_t, diff_t> Mass;

			const Chains<rank_t, diff_t>& C, D, E;
			gen_t<rank_t, diff_t> resgenC, resgenD, resgenE, boundaryCD, boundaryDE;

			ChainsBox<rank_t, diff_t> BoxCD, BoxDE;
			JunctionBox<rank_t, diff_t> Box;

			IDGeneratorCompute<rank_t, diff_t> ID;

			std::vector<rank_t> C_DE_detailedrank; //we only compute the CD_E box product and the C_DE_detailedrank

			const int level, degreeC, degreeD, degreeE, selectC, selectD, selectE;

			Eigen::PermutationMatrix<-1, -1, int> BoxDE_to_BoxCD;

			gen_t<rank_t, diff_t> box_boundary(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const ChainsBox<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);
			gen_t<rank_t, diff_t> product_bottom(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const JunctionBox<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);
			gen_t<rank_t, diff_t> product_bottom(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const std::vector<rank_t>&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);

			MasseyCompute(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int, int, int, int, int, int, int);

			void getGens();
			void box();
			void compute();
			void indeterminacy();

			///Stores Massey products and their indeterminacy
			template<typename s_rank_t, typename s_diff_t>
			friend class Mackey::Massey;
		};

		template<typename rank_t, typename diff_t>
		MasseyCompute<rank_t, diff_t>::MasseyCompute(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, const Chains<rank_t, diff_t>& E, int level, int degreeC, int degreeD, int degreeE, int selectC, int selectD, int selectE)
			: C(C), D(D), E(E), level(level), degreeC(degreeC), degreeD(degreeD), degreeE(degreeE), selectC(selectC), selectD(selectD), selectE(selectE), \
			BoxDE_to_BoxCD(box_permutation(C.rank, D.rank, E.rank, degreeC + degreeD + degreeE + 1))
		{
			getGens();
			if (Mass.isZero)
				return;
			box();
			if (!Mass.exists)
				return;
			compute();
			indeterminacy();
		}

		template<typename rank_t, typename diff_t>
		gen_t<rank_t, diff_t> MasseyCompute<rank_t, diff_t>::box_boundary(const Chains<rank_t, diff_t>& C1, const Chains<rank_t, diff_t>& C2, const ChainsBox<rank_t, diff_t>& Box, const gen_t<rank_t, diff_t>& r1, const gen_t<rank_t, diff_t>& r2, int degree1, int degree2) {
			JunctionBox<rank_t, diff_t> J_prod(Box, degree1 + degree2);
			auto prod_bottom = product_bottom(C1, C2, J_prod, r1, r2, degree1, degree2);
			IDGeneratorCompute<rank_t, diff_t> ID(level, J_prod, 1);
			auto prod_level = invRes(prod_bottom, J_prod.rank, ID.rank_level);
			auto var = ID.H_level.boundary(prod_level);
			gen_t<rank_t,diff_t> var_bottom;
			if (var.size() != 0)
				var_bottom = restriction(var, transfer(Box.rank[degree1 + degree2 + 1], level), Box.rank[degree1 + degree2 + 1]);
			return var_bottom;
		}


		template<typename rank_t, typename diff_t>
		gen_t<rank_t, diff_t>  MasseyCompute<rank_t, diff_t>::product_bottom(const Chains<rank_t, diff_t>& C1, const Chains<rank_t, diff_t>& C2, const JunctionBox<rank_t, diff_t>& Boxed, const gen_t<rank_t, diff_t>& r1, const gen_t<rank_t, diff_t>& r2, int degree1, int degree2) {
			return product_bottom(C1, C2, Boxed.detailedrank, r1, r2, degree1, degree2);
		}

		template<typename rank_t, typename diff_t>
		gen_t<rank_t, diff_t>  MasseyCompute<rank_t, diff_t>::product_bottom(const Chains<rank_t, diff_t>& C1, const Chains<rank_t, diff_t>& C2, const std::vector<rank_t>& detailedrank, const gen_t<rank_t, diff_t>& r1, const gen_t<rank_t, diff_t>& r2, int degree1, int degree2) {
			ProductGen<rank_t, diff_t> Restricted(C1, C2, degree1, degree2);
			Restricted.pad(detailedrank);
			Restricted.multiply(r1, r2);
			return Restricted.product;
		}

		template<typename rank_t, typename diff_t>
		void MasseyCompute<rank_t, diff_t>::getGens() {
			resgenC = std::move(ChainsLevelGen<rank_t, diff_t>(C, level, degreeC, selectC).res_gen);
			if (resgenC.size() == 0) {
				Mass.isZero = 1;
				return;
			}
			resgenD = std::move(ChainsLevelGen<rank_t, diff_t>(D, level, degreeD, selectD).res_gen);
			if (resgenD.size() == 0) {
				Mass.isZero = 1;
				return;
			}
			resgenE = std::move(ChainsLevelGen<rank_t, diff_t>(E, level, degreeE, selectE).res_gen);
			if (resgenE.size() == 0) {
				Mass.isZero = 1;
				return;
			}
			Mass.isZero = 0;
		}


		template<typename rank_t, typename diff_t>
		void MasseyCompute<rank_t, diff_t>::box() {
			BoxCD = ChainsBox<rank_t, diff_t>(C, D, std::min(C.maxindex + D.maxindex, degreeC + degreeD + degreeE + 2));
			boundaryCD = box_boundary(C, D, BoxCD, resgenC, resgenD, degreeC, degreeD);
			if (boundaryCD.size() == 0)
				return;
			BoxDE = ChainsBox<rank_t, diff_t>(D, E, std::min(D.maxindex + E.maxindex, degreeC + degreeD + degreeE + 2));
			boundaryDE = box_boundary(D, E, BoxDE, resgenD, resgenE, degreeD, degreeE);
			if (boundaryDE.size() == 0)
				return;
			Mass.exists = 1;
			Box = JunctionBox<rank_t, diff_t>(BoxCD, E, degreeC + degreeD + degreeE + 1);
		}

		template<typename rank_t, typename diff_t>
		void MasseyCompute<rank_t, diff_t>::compute() {
			auto factor1 = product_bottom(BoxCD, E, Box, boundaryCD, resgenE, degreeC + degreeD + 1, degreeE);

			C_DE_detailedrank = rankBox(C, BoxDE, degreeC + degreeD + degreeE + 1).second;
			gen_t<rank_t, diff_t> factor2 = BoxDE_to_BoxCD * product_bottom(C, BoxDE, C_DE_detailedrank, resgenC, boundaryDE, degreeC, degreeD + degreeE + 1);

			typename diff_t::Scalar sign = (1 - 2 * ((degreeD + degreeE) % 2));
			gen_t<rank_t, diff_t> Masseyproduct_bottom = factor1 + sign * factor2;

			ID=IDGeneratorCompute<rank_t, diff_t>(level, Box);
			Mass.boxID = std::move(ID.ID);
			Mass.isZero = ID.H_level.isZero;
			auto Masseyproduct = invRes(Masseyproduct_bottom, Box.rank, ID.rank_level);
			Mass.basis = ID.H_level.basis(Masseyproduct);
			Mass.isZero = Mass.basis.isZero();
			Mass.normalBasis = Mass.basis;
			normalize(Mass.normalBasis, ID.H_level.Groups);
			Mass.Groups = ID.H_level.Groups;
		}

		template<typename rank_t, typename diff_t>
		void MasseyCompute<rank_t, diff_t>::indeterminacy() {

			std::vector<typename rank_t::Scalar> indeterminacy_left, indeterminacy_right;

			Junction<rank_t, diff_t> J_DE(BoxDE, degreeD + degreeE + 1);
			IDGeneratorCompute<rank_t, diff_t> ID_DE(level, J_DE);

			for (int i = 0; i < ID_DE.H_level.Generators.cols(); i++) {
				gen_t<rank_t, diff_t> gen = ID_DE.H_level.Generators.col(i);
				auto resgen = restriction(gen, ID_DE.rank_level, J_DE.rank);
				gen_t<rank_t, diff_t> prod_bottom = BoxDE_to_BoxCD * product_bottom(C, BoxDE, C_DE_detailedrank, resgenC, resgen, degreeC, degreeD + degreeE + 1);
				auto prod_level = invRes(prod_bottom, Box.rank, ID.rank_level);
				auto basis=ID.H_level.basis(prod_level);
				if (basis.size() != 0 && !basis.isZero()) {
					auto o = order(basis, ID.H_level.Groups);
					if (o==0)
						indeterminacy_left.push_back(1);
					else
						indeterminacy_left.push_back(o);
				}
			}

			Junction<rank_t, diff_t> J_CD(BoxCD, degreeC + degreeD + 1);
			IDGeneratorCompute<rank_t, diff_t> ID_CD(level, J_CD);

			for (int i = 0; i < ID_CD.H_level.Generators.cols(); i++) {
				gen_t<rank_t, diff_t> gen = ID_CD.H_level.Generators.col(i);
				auto resgen = restriction(gen, ID_CD.rank_level, J_CD.rank);
				auto prod_bottom = product_bottom(BoxCD, E, Box, resgen, resgenE, degreeC+degreeD+1, degreeE);
				auto prod_level = invRes(prod_bottom, Box.rank, ID.rank_level);
				auto basis = ID.H_level.basis(prod_level);
				if (basis.size() != 0 && !basis.isZero()) {
					auto o = order(basis, ID.H_level.Groups);
					if (o == 0)
						indeterminacy_right.push_back(1);
					else
						indeterminacy_right.push_back(o);
				}
			}
			Mass.indeterminacy[0] = Eigen::Map<rank_t>(indeterminacy_left.data(), indeterminacy_left.size());
			Mass.indeterminacy[1] = Eigen::Map<rank_t>(indeterminacy_right.data(), indeterminacy_right.size());

		}
	}

	template<typename rank_t, typename diff_t>
	Massey<rank_t, diff_t>::Massey(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, const Chains<rank_t, diff_t>& E, int level, int degreeC, int degreeD, int degreeE, int selectC, int selectD, int selectE) {
		internal::MasseyCompute Comp(C, D, E, level, degreeC, degreeD, degreeE, selectC, selectD, selectE);
		*this = std::move(Comp.Mass); 
		noIndeterminacy= (indeterminacy[0].size()==0 && indeterminacy[1].size() == 0);
	}

	template<typename rank_t, typename diff_t>
	Massey<rank_t, diff_t>::Massey(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, const Chains<rank_t, diff_t>& E, int level, int degreeC, int degreeD, int degreeE)
		: Massey(C, D, E, level, degreeC, degreeD, degreeE, 0, 0, 0) {}

}
