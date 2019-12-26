#pragma once
#include "Homology.h"
#include "Automorphisms.h"
#include "MackeyFunctor.h"
#include "Z_n.h"

///@file
///@brief Contains the identification classes


namespace Mackey {

	///Stores detailed identification data for the generators in non cyclic homology groups using the Mackey functor structure
	template<typename rank_t>
	class IDGenerators {
	public:
		///The matrix type of the transfer and restriction matrices
		typedef Eigen::Matrix<typename rank_t::Scalar, -1, -1> matrix_t;
		rank_t group; ///<The homology group at our level
		rank_t group_lower; ///<The homology group at one level lower
		matrix_t Tr; ///<The transfer from one level lower to our level
		matrix_t Res; ///<The restriction from our level to one lower

		///Default Constructor
		IDGenerators() {};

		///Get a two-level Mackey functor out of this data
		MackeyFunctor<rank_t> getMackey() const;
	};

	template<typename rank_t>
	MackeyFunctor<rank_t> IDGenerators<rank_t>::getMackey() const {
		MackeyFunctor<rank_t> M;
		M.resize(2);
		M.Groups[0] = group_lower;
		M.Groups[1] = group;
		M.Tr[0] = Tr;
		M.Res[0] = Res;
		return M;
	}


///	Use the ID data to identify all possible candidates for an element
template<typename rank_t>
std::vector<rank_t> id_candidates(const rank_t& basis, const IDGenerators<rank_t>& Id_first, const IDGenerators<rank_t>& Id_second) {
	auto Mack_first = Id_first.getMackey();
	auto Mack_second = Id_second.getMackey();
	auto autom = Mack_first.automorphisms();
	std::vector<int> match;
	match.reserve(autom.first.size());
	for (int i = 0; i < autom.first.size(); i++) {
		if (Mack_second.compare(Mack_first.apply(autom.first[i], autom.second[i])))
			match.push_back(i);
	}
	std::vector<rank_t> candidates;
	candidates.reserve(match.size());
	for (int i = 0; i < match.size(); i++) {
		rank_t newbasis = autom.first[match[i]][1] * basis.transpose();
		rank_t normalbasis = normalize(newbasis,Id_second.group);
		bool flag = 0;
		for (auto& i : candidates) {
			if (i == normalbasis) {
				flag = 1;
				break;
			}
		}
		if (!flag)
			candidates.push_back(normalbasis);
	}
	return candidates;
	
}










	/////Use the ID data to identify all possible candidates for an element
	//template<typename rank_t>
	//std::vector<rank_t> id_candidates(const rank_t& basis, const IDGenerators<rank_t>& Id_first, const IDGenerators<rank_t>& Id_second) {

	//	IDGeneratorCompact<rank_t, diff_t> Id(basis, Id_first);
	//	std::vector<rank_t> candidates;

	//	rank_t maxindices = Id_second.group;
	//	rank_t minindices = maxindices;
	//	minindices.setZero();

	//	typedef typename diff_t::Scalar coeff;
	//	if constexpr (Mackey::is_finite_cyclic<coeff>()) {
	//		constexpr int order = coeff::order;
	//		for (int i = 0; i < maxindices.size(); i++) {
	//			if (maxindices[i] == 1)
	//				maxindices[i] = order - 1;
	//		}
	//	}
	//	else {
	//		for (int i = 0; i < maxindices.size(); i++) {
	//			if (maxindices[i] == 1)
	//				maxindices[i] = 0;
	//		}
	//		if (Id.torsionfree > 0) { //if torsionfree but with only one Z
	//			auto pos = find(Id_second.group, 1);
	//			minindices[pos] = maxindices[pos] = Id.torsionfree;
	//		}
	//		else if (!(Id.torsion > 0)) { //infinitely many candidates
	//			return candidates;
	//		}
	//	}
	//	candidates = DegreeConstruction(minindices, maxindices);

	//	std::vector<rank_t> better_candidates;
	//	better_candidates.reserve(candidates.size());
	//	for (const auto& i : candidates) {
	//		IDGeneratorCompact<rank_t, diff_t> Id_cand(i, Id_second);
	//		if (Id == Id_cand)
	//			better_candidates.push_back(i);
	//	}
	//	return better_candidates;
	//}

	///Contains classes and methods whose user interface is provided by the rest of the library
	namespace internal {
		///Computes the homology at given level and the identification data coming from the Mackey functor structure
		template<typename rank_t, typename diff_t>
		class IDGeneratorCompute {
		public:
			IDGenerators<rank_t> ID; ///<The identification data
			IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom); ///<Computes the identification data
		private:
			rank_t rank_level;
			Homology<rank_t, diff_t> H_level;
			IDGeneratorCompute() {};
			IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom, bool getQ);
			///Identifies the generators in case of non cyclic homology.
			void id(const Junction<rank_t, diff_t>&);

			///Compactifies the ID data
			template<typename s_rank_t>
			friend class IDGenerators;

			///Forms the input of the Multiplication Table
			template<typename s_rank_t, typename s_diff_t>
			friend class TableInput;

			///Computes products
			template<typename s_rank_t, typename s_diff_t>
			friend class GreenCompute;

			///Computes Massey products
			template<typename s_rank_t, typename s_diff_t>
			friend class MasseyCompute;

		};

		template<typename rank_t, typename diff_t>
		IDGeneratorCompute<rank_t, diff_t>::IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom, bool getQ) {
			auto current = transfer(bottom, level);
			rank_level = current.rank;
			H_level = Homology<rank_t, diff_t>(current, getQ);
			ID.group = H_level.Groups;
			if (ID.group.size() <= 1)
				return;
			else
				id(transfer(bottom, level - 1));
		}

		template<typename rank_t, typename diff_t>
		IDGeneratorCompute<rank_t, diff_t>::IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom) : IDGeneratorCompute(level, bottom, 0) {}

		template<typename rank_t, typename diff_t>
		void IDGeneratorCompute<rank_t, diff_t>::id(const Junction<rank_t, diff_t>& lower) {
			auto rank_lower = lower.rank;
			Homology<rank_t, diff_t>H_lower(lower);
			if (!H_lower.isZero) {
				ID.Tr = transfer(H_lower, H_level, rank_lower, rank_level);
				ID.Res = restriction(H_level, H_lower, rank_level, rank_lower);
				ID.group_lower = H_lower.Groups;
			}
		}
	}




//	//////////////////////////////////////////////////////////////////
/////Compactly stores usable identification data for a given generator in a non cyclic homology group using the Mackey functor structure
//
/////diff_t is used as a template to cleanly get the coefficient type.
//////////////////////////////////////////////////////////////////
//	template<typename rank_t, typename diff_t>
//	class IDGeneratorCompact {
//	public:
//		int torsion; ///<The order of the element, 1 if it's torsion free
//		int torsionfree; ///<If there is a single Z then this is the Z part of the generator. -1 otherwise
//		int restriction; ///<The order of the restriction
//		int isTransfer; ///<if the element is a transfer. -1 if we can't decide
//
//		///Constructs the identification data for our element given the detailed identification data for all generators.
//		IDGeneratorCompact(const rank_t&, const IDGenerators<rank_t, diff_t>&);
//
//		///Compares identifaction data
//		inline bool operator == (const IDGeneratorCompact& other) const {
//			return (torsion == other.torsion && torsionfree == other.torsionfree && restriction == other.restriction && isTransfer == other.isTransfer);
//		}
//
//		///Compares identifaction data
//		inline bool operator != (const IDGeneratorCompact& other) const {
//			return !(*this == other);
//		}
//	};




	//template<typename rank_t, typename diff_t>
	//IDGeneratorCompact<rank_t, diff_t>::IDGeneratorCompact(const rank_t& basis, const IDGenerators<rank_t, diff_t>& Id) {
	//	torsion = order(basis, Id.group);
	//	if (torsion == 1) {
	//		auto m = findcount(Id.group, 1);
	//		if (m.first == 1)
	//			torsionfree = basis[m.second];
	//		else
	//			torsionfree = -1;
	//	}
	//	if (Id.group_lower.size() != 0) {
	//		rank_t Res = Id.Res *basis.transpose();
	//		restriction = order(Res, Id.group_lower);
	//		isTransfer = 0;
	//		for (int i = 0; i < Id.Tr.cols(); i++) {
	//			if (basis == Id.Tr.col(i).transpose()) {
	//				isTransfer = 1;
	//				break;
	//			}
	//		}
	//		if (isTransfer == 0 && Id.group_lower.size() != 1)
	//			isTransfer = -1;
	//	}
	//	else {
	//		restriction = 0;
	//		isTransfer = 0;
	//	}
	//}



	/////Produce all generator candidates that partially match the ID of the given element
	//template<typename rank_t, typename diff_t>
	//std::vector<rank_t> id_candidates(const IDGeneratorCompact<rank_t, diff_t>& Id_element, const rank_t& group) {
	//	rank_t maxindices = group;
	//	rank_t minindices = maxindices;
	//	minindices.setZero();

	//	typedef typename diff_t::Scalar coeff;
	//	if constexpr (Mackey::is_finite_cyclic<coeff>()) {
	//		constexpr int order = coeff::order;
	//		for (int i = 0; i < maxindices.size(); i++) {
	//			if (maxindices[i] == 1)
	//				maxindices[i] = order - 1;
	//		}
	//		return DegreeConstruction(minindices, maxindices);
	//	}
	//	else {
	//		for (int i = 0; i < maxindices.size(); i++) {
	//			if (maxindices[i] == 1)
	//				maxindices[i] = 0;
	//		}
	//		if (Id_element.torsion > 0) { //if torsion then there are finitely many candidates
	//			return DegreeConstruction(minindices, maxindices);
	//		}
	//		if (Id_element.torsionfree > 0) { //if torsionfree but with only one Z
	//			auto pos = find(group, 1);
	//			minindices[pos] = maxindices[pos] = Id_element.torsionfree;
	//			return DegreeConstruction(minindices, maxindices);
	//		}
	//	}
	//	//if we have infinitely many candidates return empty
	//	return std::vector<rank_t>();
	//}
}
