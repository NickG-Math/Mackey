#pragma once
#include "Homology.h"
#include "Z_n.h"

///@file
///@brief Contains the identification classes


namespace Mackey {

	//////////////////////////////////////////////////////////////////
	///Stores detailed identification data for the generators in non cyclic homology groups using the Mackey functor structure

	///diff_t is used as a template to cleanly get the coefficient type.
	////////////////////////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	class IDGenerators {
	public:
		rank_t group; ///<The homology group at our level
		rank_t group_lower; ///<The homology group at one level lower
		std::vector<rank_t> Tr; ///<The transfers of the lower level generators to our level
		std::vector<rank_t> Res; ///<The restrictions of the generators to the lower level

		///Default Constructor
		IDGenerators() {};
	};

	//////////////////////////////////////////////////////////////////
	///Compactly stores usable identification data for a given generator in a non cyclic homology group using the Mackey functor structure

	///diff_t is used as a template to cleanly get the coefficient type.
	////////////////////////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	class IDGeneratorCompact {
	public:
		int torsion; ///<The order of the element, 1 if it's torsion free
		int torsionfree; ///<If there is a single Z then this is the Z part of the generator. -1 otherwise
		int restriction; ///<The order of the restriction
		int isTransfer; ///<if the element is a transfer. -1 if we can't decide

		///Constructs the identification data for our element given the detailed identification data for all generators.
		IDGeneratorCompact(const rank_t&, const IDGenerators<rank_t, diff_t>&);

		///Compares identifaction data
		inline bool operator == (const IDGeneratorCompact& other) const {
			return (torsion == other.torsion && torsionfree == other.torsionfree && restriction == other.restriction && isTransfer == other.isTransfer);
		}

		///Compares identifaction data
		inline bool operator != (const IDGeneratorCompact& other) const {
			return !(*this == other);
		}
	};

	template<typename rank_t, typename diff_t>
	IDGeneratorCompact<rank_t, diff_t>::IDGeneratorCompact(const rank_t& basis, const IDGenerators<rank_t, diff_t>& Id) {
		torsion = order(basis, Id.group);
		if (torsion == 1) {
			auto m = findcount(Id.group, 1);
			if (m.first == 1) {
				torsionfree = basis[m.second];
			}
			else {
				torsionfree = -1;
			}
		}
		rank_t Res = basis[0] * Id.Res[0];
		for (int i = 1; i < basis.size(); i++) {
			Res += basis[i] * Id.Res[i];
		}
		restriction = order(Res, Id.group_lower);
		isTransfer = 0;
		for (auto& i : Id.Tr) {
			if (basis == i) {
				isTransfer = 1;
				break;
			}
		}
		if (isTransfer == 0 && Id.group_lower.size() != 1) {
			isTransfer = -1;
		}
	}

	///Produce all generator candidates that match the ID of the given element
	template<typename rank_t, typename diff_t>
	std::vector<rank_t> produce_candidates(const IDGeneratorCompact<rank_t, diff_t>& Id_element, const IDGenerators<rank_t, diff_t>& Id_group) {
		rank_t maxindices = Id_group.group;
		rank_t minindices = maxindices;
		minindices.setZero();

		typedef typename diff_t::Scalar coeff;
		if constexpr (Mackey::is_finite_cyclic<coeff>()) {
			constexpr int order = coeff::order;
			for (int i = 0; i < maxindices.size(); i++) {
				if (maxindices[i] == 1)
					maxindices[i] = order - 1;
			}
			return DegreeConstruction(minindices, maxindices);
		}
		else {
			for (int i = 0; i < maxindices.size(); i++) {
				if (maxindices[i] == 1)
					maxindices[i] = 0;
			}
			if (Id_element.torsion > 0) { //if torsion then there are finitely many candidates
				return DegreeConstruction(minindices, maxindices);
			}
			if (Id_element.torsionfree > 0) { //if torsionfree but with only one Z
				auto pos = find(Id_group.group, 1);
				minindices[pos] = maxindices[pos] = Id_element.torsionfree;
				return DegreeConstruction(minindices, maxindices);
			}
		}
		//if we have infinitely many candidates...
		std::vector<rank_t> candidates;
		return candidates;
	}

	///Produce an element using the ID data
	template<typename rank_t, typename diff_t>
	std::vector<rank_t> identify(const rank_t& basis, const IDGenerators<rank_t, diff_t>& Id_first, const IDGenerators<rank_t, diff_t>& Id_second) {
		IDGeneratorCompact<rank_t, diff_t> Id(basis, Id_first);
		auto candidates = produce_candidates(Id, Id_second);
		if (candidates.empty()) {
			return candidates;
		}
		std::vector<rank_t> better_candidates;
		better_candidates.reserve(candidates.size());
		for (const auto& i:candidates) {
			IDGeneratorCompact<rank_t, diff_t> Id_cand(i, Id_second);
			if (Id == Id_cand)
				better_candidates.push_back(i);
		}
		return better_candidates;
	}

	///Contains classes and methods whose user interface is provided by the rest of the library
	namespace internal {
		///Computes the homology at given level and the identification data coming from the Mackey functor structure
		template<typename rank_t, typename diff_t>
		class IDGeneratorCompute {
		public:
			IDGenerators<rank_t, diff_t> ID; ///<The identification data
			IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom); ///<Computes the identification data
		private:
			rank_t rank_level;
			Homology<rank_t, diff_t> H_level;
			IDGeneratorCompute() {};
			IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom, bool getQ);
			///Identifies the generators in case of non cyclic homology.
			void id(const Junction<rank_t, diff_t>&);

			///Compactifies the ID data
			template<typename s_rank_t, typename s_diff_t>
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
			H_level = Homology<rank_t, diff_t>(current,getQ);
			ID.group = H_level.Groups;
			if (ID.group.size() <= 1) {
				return;
			}
			else {
				id(transfer(bottom, level - 1));
			}
		}

		template<typename rank_t, typename diff_t>
		IDGeneratorCompute<rank_t, diff_t>::IDGeneratorCompute(int level, const Junction<rank_t, diff_t>& bottom) : IDGeneratorCompute(level,bottom,0) {}

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

}
