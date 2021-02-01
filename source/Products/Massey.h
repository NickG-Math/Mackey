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
		Massey(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int, int, int, int, int, int, int, bool);
	};

	namespace internal
	{

		///General expression for the canonical basis for firstrank box secondrank using the provided names
		template<typename rank_t, typename T, typename S>
		std::vector<std::pair<T, S>> general_canonical(const rank_t& first, const rank_t& second, const std::vector<T>& firstnames, const std::vector<S>& secondnames) {
			std::vector<std::pair<T, S>> vec;
			vec.reserve(summation(first) * summation(second));
			auto s = first.size();
			auto t = second.size();
			for (int d = 0; d < t; d++) {
				auto sumsecond = summation(second, d);
				for (int c = 0; c < s; c++) {
					auto sumfirst = summation(first, c);
					for (int j = 0; j < std::min(first[c], second[d]); j++) {
						for (int i = 0; i < std::max(first[c], second[d]); i++) {
							auto pair = std::make_pair(firstnames[(i % first[c]) + sumfirst], secondnames[(i + j) % second[d] + sumsecond]);
							vec.push_back(pair);
						}
					}
				}
			}
			return vec;
		}


		///General expression for the canonical basis of (firstChain box secondChain)_i using the provided names
		template<typename rank_t, typename T, typename S>
		std::vector<std::tuple<T, S, int, int>> general_canonical_box(int degree, const std::vector<rank_t>& first, const std::vector<rank_t>& second, const std::vector<std::vector<T>>& firstnames, const std::vector<std::vector<S>>& secondnames) {
			std::vector<std::tuple<T, S, int, int>> vec;
			vec.reserve((summation(rankBox(first, second, degree).first)));
			auto s = first.size();
			auto t = second.size();
			long sum = 0;
			for (int i = 0; i < s; i++) {
				int j = degree - i;
				if (j >=0 && j <t){
					auto u = general_canonical(first[i], second[j], firstnames[i], secondnames[j]);
					for (const auto& pair : u)
						vec.push_back(std::make_tuple(pair.first, pair.second, i, j));
				}
			}
			return vec;
		}



		///General expression for the canonical basis of firstChain box secondChain using the provided names
		template<typename rank_t, typename T, typename S>
		std::vector<std::vector<std::tuple<T, S, int, int>>> general_canonical_box(const std::vector<rank_t>& first, const std::vector<rank_t>& second, const std::vector<std::vector<T>>& firstnames, const std::vector<std::vector<S>>& secondnames) {
			std::vector<std::vector<std::tuple<T, S, int, int>>> vec(first.size() + second.size() - 1);
			for (int i = 0; i < vec.size(); i++)
				vec[i] = general_canonical_box(i, first, second, firstnames, secondnames);
			return vec;
		}


		///Produces a nicer version of general canonical when used for triple products
		std::array<long, 6> triple_box_appender(const std::tuple<std::tuple<long, long, int, int>, long, int, int>& a) {
			auto u = std::get<0>(a);
			long first = std::get<0>(u);
			long second = std::get<1>(u);
			long firstdeg = std::get<2>(u);
			long seconddeg = std::get<3>(u);
			long third = std::get<1>(a);
			long thirddeg = std::get<3>(a); // std::get<2>(i) is firstseconddeg=firstdeg+seconddeg hence need not be recorded
			return { first,second,third,firstdeg,seconddeg,thirddeg };
		}

		///Produces a nicer version of general canonical when used for triple products
		std::array<long, 6> triple_box_appender(const std::tuple<long, std::tuple<long, long, int, int>, int, int>& a) {
			long first = std::get<0>(a);
			auto u = std::get<1>(a);
			long second = std::get<0>(u);
			long third = std::get<1>(u);
			long seconddeg = std::get<2>(u);
			long thirddeg = std::get<3>(u);
			long firstdeg = std::get<2>(a);
			return { first,second,third,firstdeg,seconddeg,thirddeg };
		}

		///Forms the permutation from basis (C box D) box E to C box (D box E) using names provided by general_canonical
		template<typename T, typename S>
		decltype(auto) triple_box_permutation(const std::vector<T>& a, const std::vector<S>& b) { 
			std::vector<std::array<long, 6>> fixed_a, fixed_b;
			fixed_a.reserve(a.size());
			fixed_b.reserve(b.size());
			for (const auto& i : a)
				fixed_a.push_back(triple_box_appender(i));
			for (const auto& i : b)
				fixed_b.push_back(triple_box_appender(i));
			return changebasis<long>(fixed_a, fixed_b);
		}


		///Returns the names for an equivariant basis
		template<typename rank_t>
		std::vector<long> equivariant_base_names(const rank_t& first) {
			std::vector<long> vec(summation(first));
			std::iota(vec.begin(), vec.end(), 0);
			return vec;
		}

		///Forms the permutation allowing us to go from C box (D box E) to (C box D) box E.
		template<typename rank_t>
		decltype(auto) box_permutation(const std::vector<rank_t>& C, const std::vector<rank_t>& D, const std::vector<rank_t>& E, int degree) {
			std::vector<std::vector<long>> Cnames, Dnames, Enames;
			Cnames.reserve(C.size()); Dnames.reserve(D.size()); Enames.reserve(E.size());

			for (const auto& i : C)
				Cnames.push_back(equivariant_base_names(i));
			for (const auto& i : D)
				Dnames.push_back(equivariant_base_names(i));
			for (const auto& i : E)
				Enames.push_back(equivariant_base_names(i));

			auto CD_names = general_canonical_box(C, D, Cnames, Dnames);
			auto CD = rankBox(C, D);
			auto CD_E = general_canonical_box(degree, CD, E, CD_names, Enames);

			auto DE_names = general_canonical_box(D, E, Dnames, Enames);
			auto DE = rankBox(D, E);
			auto C_DE = general_canonical_box(degree, C, DE, Cnames, DE_names);
			
			auto p = triple_box_permutation(CD_E, C_DE);
			if (p.empty())
				p = { 0 };
			return Eigen::PermutationMatrix<-1, -1, long>(Eigen::Map<Eigen::Matrix<long, -1, 1>>(p.data(), p.size()));
		}


		///Computes Massey products and their indeterminacy
		template<typename rank_t, typename diff_t>
		struct MasseyCompute {
			Massey<rank_t, diff_t> Mass;

			const Chains<rank_t, diff_t>& C, D, E;
			gen_t<rank_t, diff_t> resgenC, resgenD, resgenE, boundaryCD, boundaryDE;

			ChainsBox<rank_t, diff_t> BoxCD, BoxDE;
			JunctionBox<rank_t, diff_t> Box;

			IDGeneratorCompute<rank_t, diff_t> ID;

			std::vector<long> C_DE_detailedrank; //we only compute the CD_E box product and the C_DE_detailedrank

			const int level, degreeC, degreeD, degreeE, selectC, selectD, selectE;
			Eigen::PermutationMatrix<-1, -1, long> BoxDE_to_BoxCD;

			gen_t<rank_t, diff_t> box_boundary(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const ChainsBox<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);
			gen_t<rank_t, diff_t> product_bottom(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const JunctionBox<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);
			gen_t<rank_t, diff_t> product_bottom(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const std::vector<long>&, const gen_t<rank_t, diff_t>&, const gen_t<rank_t, diff_t>&, int, int);

			MasseyCompute(const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, const Chains<rank_t, diff_t>&, int, int, int, int, int, int, int, bool);

			void getGens();
			void box();
			void compute();
			void indeterminacy();

			///Stores Massey products and their indeterminacy
			template<typename, typename>
			friend class Mackey::Massey;
		};

		template<typename rank_t, typename diff_t>
		MasseyCompute<rank_t, diff_t>::MasseyCompute(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, const Chains<rank_t, diff_t>& E, int level, int degreeC, int degreeD, int degreeE, int selectC, int selectD, int selectE, bool do_indeter)
			: C(C), D(D), E(E), level(level), degreeC(degreeC), degreeD(degreeD), degreeE(degreeE), selectC(selectC), selectD(selectD), selectE(selectE)
		{
			BoxDE_to_BoxCD = box_permutation(C.rank, D.rank, E.rank, degreeC + degreeD + degreeE + 1);
			getGens();
			if (Mass.isZero)
				return;
			box();
			if (!Mass.exists)
				return;
			compute(); //needs to happen before indeterminacy as it computes the homology in the triple box
			if (do_indeter)
				indeterminacy();
		}


		template<typename rank_t, typename diff_t>
		gen_t<rank_t, diff_t> MasseyCompute<rank_t, diff_t>::box_boundary(const Chains<rank_t, diff_t>& C1, const Chains<rank_t, diff_t>& C2, const ChainsBox<rank_t, diff_t>& Box, const gen_t<rank_t, diff_t>& r1, const gen_t<rank_t, diff_t>& r2, int degree1, int degree2) {
			JunctionBox<rank_t, diff_t> J_prod(Box, degree1 + degree2);
			auto prod_bottom = product_bottom(C1, C2, J_prod, r1, r2, degree1, degree2);
			IDGeneratorCompute<rank_t, diff_t> ID(level, J_prod, 1);
			auto prod_level = invRes(prod_bottom, J_prod.rank, ID.rank_level);
			auto var = ID.H_level.boundary(prod_level);
			gen_t<rank_t, diff_t> var_bottom;

			if (var.size() != 0)
				var_bottom = restriction(var, transfer(Box.rank[degree1 + degree2 + 1], level), Box.rank[degree1 + degree2 + 1]);
			return var_bottom;
		}


		template<typename rank_t, typename diff_t>
		gen_t<rank_t, diff_t>  MasseyCompute<rank_t, diff_t>::product_bottom(const Chains<rank_t, diff_t>& C1, const Chains<rank_t, diff_t>& C2, const JunctionBox<rank_t, diff_t>& Boxed, const gen_t<rank_t, diff_t>& r1, const gen_t<rank_t, diff_t>& r2, int degree1, int degree2) {
			return product_bottom(C1, C2, Boxed.detailedrank, r1, r2, degree1, degree2);
		}

		template<typename rank_t, typename diff_t>
		gen_t<rank_t, diff_t>  MasseyCompute<rank_t, diff_t>::product_bottom(const Chains<rank_t, diff_t>& C1, const Chains<rank_t, diff_t>& C2, const std::vector<long>& detailedrank, const gen_t<rank_t, diff_t>& r1, const gen_t<rank_t, diff_t>& r2, int degree1, int degree2) {
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
			C_DE_detailedrank = rankBox(C.rank, BoxDE.rank, degreeC + degreeD + degreeE + 1).second;
			Box = JunctionBox<rank_t, diff_t>(BoxCD, E, degreeC + degreeD + degreeE + 1);
		}

		template<typename rank_t, typename diff_t>
		void MasseyCompute<rank_t, diff_t>::compute() {
			auto factor1 = product_bottom(BoxCD, E, Box, boundaryCD, resgenE, degreeC + degreeD + 1, degreeE);
			gen_t<rank_t, diff_t> factor2 = BoxDE_to_BoxCD * product_bottom(C, BoxDE, C_DE_detailedrank, resgenC, boundaryDE, degreeC, degreeD + degreeE + 1);

			typename diff_t::Scalar sign;
			if (degreeC % 2)
				sign = 1;
			else
				sign = -1;
			gen_t<rank_t, diff_t> Masseyproduct_bottom = factor1 + sign * factor2;

			ID = IDGeneratorCompute<rank_t, diff_t>(level, Box);
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
				auto basis = ID.H_level.basis(prod_level);
				if (basis.size() != 0 && !basis.isZero()) {
					auto o = order(basis, ID.H_level.Groups);
					if (o == 0)
						indeterminacy_left.push_back(1);
					else
						indeterminacy_left.push_back(o);
				}
			}
			Mass.noIndeterminacy = indeterminacy_left.empty();
			Junction<rank_t, diff_t> J_CD(BoxCD, degreeC + degreeD + 1);
			IDGeneratorCompute<rank_t, diff_t> ID_CD(level, J_CD);

			for (int i = 0; i < ID_CD.H_level.Generators.cols(); i++) {
				gen_t<rank_t, diff_t> gen = ID_CD.H_level.Generators.col(i);
				auto resgen = restriction(gen, ID_CD.rank_level, J_CD.rank);
				auto prod_bottom = product_bottom(BoxCD, E, Box, resgen, resgenE, degreeC + degreeD + 1, degreeE);
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

			Mass.noIndeterminacy = Mass.noIndeterminacy && indeterminacy_right.empty();
			Mass.indeterminacy[0] = Eigen::Map<rank_t>(indeterminacy_left.data(), indeterminacy_left.size());
			Mass.indeterminacy[1] = Eigen::Map<rank_t>(indeterminacy_right.data(), indeterminacy_right.size());
		}
	}

	template<typename rank_t, typename diff_t>
	Massey<rank_t, diff_t>::Massey(const Chains<rank_t, diff_t>& C, const Chains<rank_t, diff_t>& D, const Chains<rank_t, diff_t>& E, int level, int degreeC, int degreeD, int degreeE, int selectC, int selectD, int selectE, bool do_indeter) {
		internal::MasseyCompute Comp(C, D, E, level, degreeC, degreeD, degreeE, selectC, selectD, selectE, do_indeter);
		*this = std::move(Comp.Mass);
	}
}
