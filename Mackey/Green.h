#pragma once
#include "Levels.h"
#include "Box.h"
#include "Homology.h"
#include <numeric> //for gcd

///@file
///@brief Contains the class and methods for multiplying generators.

namespace Mackey {

	/// The result of multiplying two generators in a Green functor
	template<typename rank_t, typename diff_t>
	class Green {
	public:
		
		rank_t basis;///<Expresses the product of two generators as a linear combination of the generators of the box product.
		rank_t normalBasis;///<Same as basis, but normalized to not have signs.
		rank_t Groups;///<The groups of the tensor product.
		bool isZero;///<1 if the result is 0

		///Computes the generators genC, genD of two Chains C,D using their degrees and selections, boxes the C,D and computes genC*genD in the box product.
		Green(Chains<rank_t, diff_t>&, Chains<rank_t, diff_t>&, int, int, int, int, int);

		///Same as the other constructor but with default selections 0,0.
		Green(Chains<rank_t, diff_t>&, Chains<rank_t, diff_t>&, int, int, int);

	private:
		Chains<rank_t, diff_t>& C, D;
		rank_t rankC_level, rankD_level, rankBoxed_bottom;

		typedef Eigen::Matrix<typename diff_t::Scalar, -1, 1> gen_t;
		gen_t genC, genD;

		const int level, degreeC, degreeD;
		int padleft, padright;
		Junction<rank_t, diff_t> Boxed_Level;

		/// Computes the generators from their degrees and the selections
		void getGens(int, int);
		/// Computes the tensor product of Chains C,D
		void box();
		/// Computes the product of generators genC, genD by combining the results of getGens and box
		void compute();


		///Rechooses the generators with new selections and computes their product. Not currently used
		void rechoose(int selectC, int selectD) {
			getGens(selectC, selectD);
			compute();
		}

	};


	template<typename rank_t, typename diff_t>
	Green<rank_t, diff_t> ::Green(Chains<rank_t, diff_t>& C, Chains<rank_t, diff_t>& D, int level, int degreeC, int degreeD, int selectC, int selectD)
		: C(C), D(D), level(level), degreeC(degreeC), degreeD(degreeD)
	{
		getGens(selectC, selectD);
		if (isZero) {
			return;
		}
		box();
		compute();
	}

	template<typename rank_t, typename diff_t>
	Green<rank_t, diff_t>::Green(Chains<rank_t, diff_t>& C, Chains<rank_t, diff_t>& D, int level, int degreeC, int degreeD) :Green(C, D, level, degreeC, degreeD, 0, 0) {}
	   
	template<typename rank_t, typename diff_t>
	void Green<rank_t, diff_t>::getGens(int selectC, int selectD) {
		Junction<rank_t, diff_t> JC(C, degreeC);
		Junction<rank_t, diff_t> JC_Level = transfer(JC, level);
		Homology<rank_t, diff_t> HC(JC_Level);
		if (HC.Groups.size() > 1) {
			genC = HC.Generators.col(selectC);
		}
		else {
			genC = HC.Generators.col(0);
		}
		rankC_level = std::move(JC_Level.rank);

		Junction<rank_t, diff_t> JD(D, degreeD);
		Junction<rank_t, diff_t> JD_Level = transfer(JD, level);
		Homology<rank_t, diff_t> HD(JD_Level);

		if (HD.Groups.size() > 1) {
			genD = HD.Generators.col(selectD);
		}
		else {
			genD = HD.Generators.col(0);
		}
		rankD_level = std::move(JD_Level.rank);

		isZero = HC.isZero & HD.isZero;
	}

	template<typename rank_t, typename diff_t>
	void Green<rank_t, diff_t>::box() {
		JunctionBox<rank_t, diff_t> Boxed(C, D, degreeC + degreeD);
		rankBoxed_bottom = Boxed.rank;
		Boxed_Level = transfer(Boxed, level);
		padleft = padright = 0;
		for (int i = 0; i <= degreeC - 1; i++) {
			padleft += summation(Boxed.detailedrank[i]);
		}
		for (int i = degreeC + 1; i <= degreeC + degreeD; i++) {
			padright += summation(Boxed.detailedrank[i]);
		}
	}

	template<typename rank_t, typename diff_t>
	void Green<rank_t, diff_t>::compute() {

		gen_t resGenC = restriction(genC, rankC_level, C.rank[degreeC]);
		gen_t resGenD = restriction(genD, rankD_level, D.rank[degreeD]);

		gen_t leftConvProduct(resGenC.size() * resGenD.size());
		for (int i = 0; i < resGenD.size(); i++) {
			leftConvProduct.segment(i * resGenC.size(), resGenC.size()) = resGenC * resGenD(i);
		}
		ChangeBasis<rank_t> permute(C.rank[degreeC], D.rank[degreeD]);
		gen_t canonProduct = permute.LefttoCanon.inverse() * leftConvProduct;
		gen_t restrictedProduct(padleft + canonProduct.size() + padright);
		restrictedProduct.setZero();
		restrictedProduct.segment(padleft, canonProduct.size()) = canonProduct;

		gen_t product = invRes(restrictedProduct, rankBoxed_bottom, Boxed_Level.rank);
		Homology<rank_t, diff_t> H(Boxed_Level);

		basis = H.basis(product);
		Groups = std::move(H.Groups);
		normalBasis.resize(basis.size());
		for (int i = 0; i < basis.size(); i++) {
			if (Groups(i) != 1 && basis(i) != 0) {
				normalBasis(i) = std::gcd(basis(i), Groups(i));
			}
			else {
				normalBasis(i) = abs(basis(i));
			}
		}
	}
}
