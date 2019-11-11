#pragma once
#include "Wrapper.h"
#include "Chains.h" 
#include "General.h" 

///@file
///@brief Contains the methods for transfering, restricting and Weyl-group-acting.

namespace Mackey {

	///Transfer the rank to given level.
	template<typename Scalar>
	Eigen::Matrix<Scalar, 1, -1> transfer(const Eigen::Matrix<Scalar, 1, -1>& rank, int level) {
		Eigen::Matrix<Scalar, 1, -1> transferred = rank;
		auto cutoff = intexp(power - level);
		for (int i = 0; i < transferred.size(); i++) //when Eigen gets iterators we can replace this with auto
		{
			if (transferred(i) > cutoff)
				transferred(i) = cutoff;
		}
		return transferred;
	}

	///Transfer the differential to given level. We need the ranks both at the original level and the given level, for both domain and range.
	template<typename rank_t, typename diff_t>
	diff_t transfer(const diff_t& diff, const rank_t& domain, rank_t& domain_top, const rank_t& range, rank_t& range_top, int level) {
		auto n = intexp(power - level);
		if (diff.size() == 0) {
			return diff;
		}
		diff_t transfer(summation(range_top), summation(domain_top));
		std::vector<int> keep;
		keep.reserve(range.size() * n);
		int sum = 0;
		int j = 0;
		for (int i = 0; i < range.size(); i++) {
			while (j < n + sum + range(i) && j < diff.rows()) {
				if (!(n + sum <= j && j <= range(i) + sum - 1)) {
					keep.push_back(j);
				}
				j++;
			}
			sum += range(i);
		}
		diff_t reduceddiff = KeepRow(diff, keep);
		int track = 0;
		int tracktransfer = 0;
		for (int i = 0; i < domain.size(); i++)
		{
			auto limit = domain(i) + track;
			if (limit - n < track + 1) {
				for (int j = track; j < limit; j++) {
					transfer.col(tracktransfer) = reduceddiff.col(j);
					tracktransfer++;
				}
			}
			else {
				for (int j = track; j < std::min(track + n, limit - n); j++) {
					Eigen::Matrix<typename diff_t::Scalar, -1, 1> columnsum = reduceddiff.col(j);
					for (int k = j + n; k < limit; k += n) {
						columnsum += reduceddiff.col(k);
					}
					transfer.col(tracktransfer) = columnsum;
					tracktransfer++;
				}
			}
			track += domain(i);
		}
		return transfer;

	}
	/////////////////////////////////////////////////
	/// Storage of various levels of Junction/Chains

	/// T is any type with a transfer function. Example: Junction, Chains.
	//////////////////////////////////////
	template<typename T>
	class Levels {
	public:
		std::vector<T> level;
		Levels(T& bottom) {
			level.resize(power + 1);
			level[0] = std::move(bottom);
			for (int i = 1; i < power + 1; i++) {
				level[i] = transfer(level[0], i);
			}
		}
	};

	///Transfer Junction to the desired level.
	template<typename rank_t, typename diff_t>
	Junction<rank_t, diff_t> transfer(const Junction<rank_t, diff_t>& J, int level) {
		rank_t rankIn_level = transfer(J.rankIn, level);
		rank_t rank_level = transfer(J.rank, level);
		rank_t rankOut_level = transfer(J.rankOut, level);
		diff_t diffIn_level = transfer(J.diffIn, J.rankIn, rankIn_level, J.rank, rank_level, level);
		diff_t diffOut_level = transfer(J.diffOut, J.rank, rank_level, J.rankOut, rankOut_level, level);
		Junction<rank_t, diff_t> J_Level(rank_level, rankOut_level, rankIn_level, diffOut_level, diffIn_level);
		return J_Level;
	}

	///Transfer Chains to the desired level.
	template<typename rank_t, typename diff_t>
	Chains<rank_t, diff_t> transfer(const Chains<rank_t, diff_t>& C, int level) {
		Chains<rank_t, diff_t> C_Level(C.rank.size());
		for (const auto& i : C.rank) {
			C_Level.rank.push_back(transfer(i, level));
		}
		C_Level.diff.resize(C.diff.size());
		for (int i = 1; i <= C.maxindex; i++) {
			C_Level.diff[i] = transfer(C.diff[i], C.rank[i], C_Level.rank[i], C.rank[i - 1], C_Level.rank[i - 1], level);
		}
		return C_Level;
	}


	///Transfer generator to level given the ranks at the original level (domain) and the target level (range).
	template<typename rank_t, typename Derived>
	Derived transfer(const Eigen::MatrixBase<Derived>& generator, const rank_t& domain, const rank_t& range) {
		Derived transferred(summation(range));
		int trackdom = 0;
		int trackran = 0;
		for (int i = 0; i < domain.size(); i++)
		{
			if (domain(i) == range(i)) {
				transferred.segment(trackran, range(i)) = prime * generator.segment(trackdom, range(i));
			}
			else {
				transferred.segment(trackran, range(i)) = generator.segment(trackdom, range(i));
				for (int k = 1; k < domain(i) / range(i); k++) {
					transferred.segment(trackran, range(i)) += generator.segment(trackdom + k * range(i), range(i));
				}
			}
			trackdom += domain(i);
			trackran += range(i);
		}
		return transferred;
	}


	///Restrict generator to level given the ranks at the original level (domain) and the target level (range).
	template<typename rank_t, typename gen_t>
	gen_t restriction(const gen_t& generator, const rank_t& domain, const rank_t& range)
	{
		gen_t restricted(summation(range));
		int trackdom = 0;
		int trackran = 0;
		for (int i = 0; i < domain.size(); i++) {
			for (int j = 0; j < range(i) / domain(i); j++) {
				restricted.segment(trackran + j * domain(i), domain(i)) = generator.segment(trackdom, domain(i));
			}
			trackdom += domain(i);
			trackran += range(i);
		}
		return restricted;
	}

	
	/////////////////////////////////////////////
	///A variant of \ref rotate "rotate" for segments of Eigen vectors
	///
	///The reason we do shiftsegment(V,a,b) as opposed to rotate(V.segment(a,b)) is explained in "In which cases do functions taking a plain Matrix or Array argument fail?" of Eigen documentation. The trick there results in segmentation faults...
	//////////////////////////////////////////////
	template<typename Derived, typename T>
	void shiftsegment(Eigen::MatrixBase<Derived>& V, int tracker, T rankvalue)
	{
		for (int i = tracker; i < tracker + rankvalue - 1; i++) {
			std::swap(V(i), V(i + 1));
		}
	}

	///Compute the Weyl group action on a generator given the rank of the group it lives in.
	template<typename rank_t, typename gen_t>
	gen_t action(const gen_t& generator, const rank_t& rank)
	{
		gen_t Weyl = generator;
		int tracker = 0;
		for (int i = 0; i < rank.size(); i++) {
			shiftsegment(Weyl, tracker, rank(i));
			tracker += rank(i);
		}
		return Weyl;
	}

	///The inverse of the restriction function on a generator given the ranks at the original level (domain) and the target level (range).
	///
	/// Recall that free Mackey functors have injective restrictions, so we only need the generator to be in the image.
	template<typename rank_t, typename gen_t>
	gen_t invRes(const gen_t& generator, const rank_t& domain, const rank_t& range)
	{
		gen_t unrestricted(summation(range));
		int trackdom = 0, trackran = 0;
		for (int i = 0; i < range.size(); i++) {
			unrestricted.segment(trackran, range(i)) = generator.segment(trackdom, range(i));
			trackdom += domain(i);
			trackran += range(i);
		}
		return unrestricted;
	}
}
