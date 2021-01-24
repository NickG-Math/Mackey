#pragma once
#include "Wrapper.h"
#include "Homology.h" 
#include "General.h" 

///@file
///@brief Contains the methods for transfering, restricting and Weyl-group-acting.

namespace Mackey {

	///Transfer the rank to given level.
	template<typename Scalar>
	Eigen::Matrix<Scalar, 1, -1> transfer(const Eigen::Matrix<Scalar, 1, -1>& rank, int level) {
		Eigen::Matrix<Scalar, 1, -1> transferred(rank.size());
		Scalar cutoff = intexp(power - level);
		for (int i = 0; i < transferred.size(); i++) //when Eigen gets iterators we can replace this with auto
			transferred[i] = std::min(rank[i], cutoff);
		return transferred;
	}

	///Transfer the differential to given level. We need the ranks both at the original level and the given level, for both domain and range.
	template<typename rank_t, typename diff_t>
	diff_t transfer(const diff_t& diff, const rank_t& domain, rank_t& domain_top, const rank_t& range, rank_t& range_top, int level) {
		if (diff.size() == 0)
			return diff;

		auto n = intexp(power - level);

		typedef typename std::conditional<SFINAE::is_Sparse<diff_t>::value, spm_t<diff_t>, diff_t>::type coltype;
		coltype reduceddiff;
		typename diff_t::StorageIndex track = 0;
		typename diff_t::StorageIndex tracktransfer = 0;
		coltype transfer(summation(range_top), summation(domain_top));

		if constexpr (SFINAE::is_Sparse<diff_t>::value) {
			std::vector<typename diff_t::StorageIndex> keep(diff.rows(), -1);
			typename diff_t::StorageIndex sum = 0;
			typename diff_t::StorageIndex j = 0;
			long counter = 0;
			for (int i = 0; i < range.size(); i++) {
				while (j < n + sum + range(i) && j < diff.rows())
				{
					if (!(n + sum <= j && j <= range(i) + sum - 1)) {
						keep[j] = counter;
						counter++;
					}
					j++;
				}
				sum += range(i);
			}
			auto u = keep_row_triplets(diff, keep);
			reduceddiff.resize(counter, diff.cols());
			reduceddiff.setFromTriplets(u.begin(), u.end());
		}
		else {
			std::vector<typename diff_t::StorageIndex> keep;
			keep.reserve(range.size() * n);
			typename diff_t::StorageIndex sum = 0;
			typename diff_t::StorageIndex j = 0;
			for (int i = 0; i < range.size(); i++) {
				while (j < n + sum + range(i) && j < diff.rows()) {
					if (!(n + sum <= j && j <= range(i) + sum - 1))
						keep.push_back(j);
					j++;
				}
				sum += range(i);
			}
			reduceddiff = KeepRow(diff, keep);
		}

		if constexpr (SFINAE::is_Sparse<diff_t>::value) { //estimate the nonzero elements in each column
			std::vector<long> cols;
			cols.reserve(diff.cols());
			for (int i = 0; i < domain.size(); i++)
			{
				auto limit = domain(i) + track;
				if (limit - n < track + 1) {
					for (auto j = track; j < limit; j++)
						cols.push_back(reduceddiff.col(j).nonZeros());
				}
				else {
					for (auto j = track; j < std::min(track + n, limit - n); j++) {
						auto counter = reduceddiff.col(j).nonZeros();
						for (auto k = j + n; k < limit; k += n)
							counter += reduceddiff.col(j).nonZeros();
						cols.push_back(counter);
					}
				}
				track += domain(i);
			}
			track = 0;
			transfer.reserve(cols);
		}

		for (int i = 0; i < domain.size(); i++)
		{
			auto limit = domain(i) + track;
			if (limit - n < track + 1) {
				for (auto j = track; j < limit; j++) {
					transfer.col(tracktransfer) = reduceddiff.col(j);
					tracktransfer++;
				}
			}
			else {
				for (auto j = track; j < std::min(track + n, limit - n); j++) {
					transfer.col(tracktransfer) = reduceddiff.col(j);
					for (auto k = j + n; k < limit; k += n)
						transfer.col(tracktransfer) += reduceddiff.col(k);
					tracktransfer++;
				}
			}
			track += domain(i);
		}
		if constexpr (SFINAE::is_Sparse<diff_t>::value) 
			transfer.prune(0,0); //remove 0 elements
		return transfer;
	}

	/////////////////////////////////////////////////
	/// Storage of various levels of Junction/Chains

	/// T is any type with a transfer function. Example: Junction, Chains.
	//////////////////////////////////////
	template<typename T>
	class Levels {
	public:
		std::vector<T> level; ///<The various levels.

		///Transfer the bottom and get all levels.
		Levels(T& bottom) {
			level.resize(power + 1);
			level[0] = std::move(bottom);
			for (int i = 1; i < power + 1; i++)
				level[i] = transfer(level[0], i);
		}
	};

	///Transfer Junction to the desired level.
	template<typename rank_t, typename diff_t>
	Junction<rank_t, diff_t> transfer(const Junction<rank_t, diff_t>& J, int level) {
		auto rankIn_level = transfer(J.rankIn, level);
		auto rank_level = transfer(J.rank, level);
		auto rankOut_level = transfer(J.rankOut, level);
		auto diffIn_level = transfer(J.diffIn, J.rankIn, rankIn_level, J.rank, rank_level, level);
		auto diffOut_level = transfer(J.diffOut, J.rank, rank_level, J.rankOut, rankOut_level, level);
		Junction<rank_t, diff_t> J_Level(rank_level, rankOut_level, rankIn_level, diffOut_level, diffIn_level);
		return J_Level;
	}

	///Transfer Chains to the desired level.
	template<typename rank_t, typename diff_t>
	Chains<rank_t, diff_t> transfer(const Chains<rank_t, diff_t>& C, int level) {
		Chains<rank_t, diff_t> C_Level(C.rank.size());
		for (const auto& i : C.rank)
			C_Level.rank.push_back(transfer(i, level));
		C_Level.diff.resize(C.diff.size());
		for (int i = 1; i <= C.maxindex; i++)
			C_Level.diff[i] = transfer(C.diff[i], C.rank[i], C_Level.rank[i], C.rank[i - 1], C_Level.rank[i - 1], level);
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
			if (domain(i) == range(i))
				transferred.segment(trackran, range(i)) = static_cast<Scalar_t<Derived>>(prime)* generator.segment(trackdom, range(i));
			else {
				transferred.segment(trackran, range(i)) = generator.segment(trackdom, range(i));
				for (int k = 1; k < domain(i) / range(i); k++)
					transferred.segment(trackran, range(i)) += generator.segment(trackdom + k * range(i), range(i));
			}
			trackdom += domain(i);
			trackran += range(i);
		}
		return transferred;
	}

	///Restrict generator to level given the ranks at the original level (domain) and the target level (range).
	template<typename rank_t, typename Derived>
	Derived restriction(const Eigen::MatrixBase<Derived>& generator, const rank_t& domain, const rank_t& range)
	{
		Derived restricted(summation(range));
		int trackdom = 0;
		int trackran = 0;
		for (int i = 0; i < domain.size(); i++) {
			for (int j = 0; j < range(i) / domain(i); j++)
				restricted.segment(trackran + j * domain(i), domain(i)) = generator.segment(trackdom, domain(i));
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
		for (int i = tracker; i < tracker + rankvalue - 1; i++)
			std::swap(V(i), V(i + 1));
	}

	///Compute the Weyl group action on a generator given the rank of the group it lives in.
	template<typename rank_t, typename Derived>
	Derived action(const Eigen::MatrixBase<Derived>& generator, const rank_t& rank)
	{
		Derived Weyl = generator;
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
	template<typename rank_t, typename T>
	T invRes(const T& generator, const rank_t& domain, const rank_t& range)
	{
		T unrestricted(summation(range));
		long trackdom = 0, trackran = 0;
		for (int i = 0; i < range.size(); i++) {
			unrestricted.segment(trackran, range(i)) = generator.segment(trackdom, range(i));
			trackdom += domain(i);
			trackran += range(i);
		}
		return unrestricted;
	}


	///Writing the transfer of each generator in terms of the generators in the image.
	template<typename rank_t, typename diff_t>
	mat_t<rank_t> transfer(const Homology<rank_t, diff_t>& low, const Homology<rank_t, diff_t>& high, const rank_t& rank_low, const rank_t& rank_high) {
		mat_t<rank_t> Tr(high.Groups.size(), low.Groups.size());
		for (long i = 0; i < low.Generators.cols(); i++) {
			gen_t<rank_t, diff_t> generator = low.Generators.col(i);
			Tr.col(i) = high.basis(transfer(generator, rank_low, rank_high)).transpose();
		}
		return Tr;
	}

	///Writing the restriction of each generator in terms of the generators in the image.
	template<typename rank_t, typename diff_t>
	mat_t<rank_t> restriction(const Homology<rank_t, diff_t>& high, const Homology<rank_t, diff_t>& low, const rank_t& rank_high, const rank_t& rank_low) {
		mat_t<rank_t> Res(low.Groups.size(), high.Groups.size());
		for (long i = 0; i < high.Generators.cols(); i++) {
			gen_t<rank_t, diff_t> generator = high.Generators.col(i);
			Res.col(i) = low.basis(restriction(generator, rank_high, rank_low)).transpose();
		}
		return Res;
	}

	///Writing the Weyl group action on each generator in terms of the other generators.
	template<typename rank_t, typename diff_t>
	mat_t<rank_t> action(const Homology<rank_t, diff_t>& H, const rank_t& rank) {
		mat_t<rank_t> Weyl(H.Groups.size(), H.Groups.size());
		for (long i = 0; i < H.Generators.cols(); i++) {
			gen_t<rank_t, diff_t> generator = H.Generators.col(i);
			Weyl.col(i) = H.basis(action(generator, rank)).transpose();
		}
		return Weyl;
	}
}
