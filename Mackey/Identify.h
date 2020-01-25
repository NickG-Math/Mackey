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
	bool operator==(const IDGenerators<rank_t>& I, const IDGenerators<rank_t>& J) {
		return I.group.size() == J.group.size() && I.group == J.group && I.group_lower.size() == J.group_lower.size() && I.group_lower == J.group_lower && I.Tr == J.Tr && I.Res == J.Res;
	}



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
		normalize(newbasis,Id_second.group); 
		bool flag = 0;
		for (auto& i : candidates) {
			if (i == newbasis) {
				flag = 1;
				break;
			}
		}
		if (!flag)
			candidates.push_back(newbasis);
	}
	return candidates;
}


///	Find if given elements can all be distinguished i.e. they have pairwise disjoint candidate sets
template<typename rank_t>
bool distinguish(const std::vector<rank_t>& elements, const IDGenerators<rank_t>& Id) {

	if (Id.group.size()==1){
		for (int i = 0; i < elements.size(); i++) {
			for (int j = i + 1; j < elements.size(); j++) {
				if (elements[i]==elements[j])
					return 0;
			}
		}
		return 1;
	}
	  
	std::vector<std::vector<rank_t>> all_candidates;
	all_candidates.reserve(elements.size());
	for (const auto& i : elements) {
		all_candidates.push_back(id_candidates(i,Id,Id));
	}
	for (int i = 0; i < all_candidates.size(); i++) {
		for (int j = i + 1; j < all_candidates.size(); j++) {
			if (find(all_candidates[i], all_candidates[j][0]) != -1)
				return 0;
		}
	}
	return 1;
}


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
			H_level = Homology<rank_t, diff_t>(current,getQ);

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

}
