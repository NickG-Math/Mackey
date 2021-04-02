#pragma once
#include "../Levels.hpp"

///@file
///@brief Contains the methods for transfering, restricting and Weyl-group-acting.



namespace mackey {

	namespace implementation_details {
		///Raise power to given exponent
		template<typename T, typename t>
		t intexp(const t exponent) {
			if constexpr (T::prime == 2)
				return 1 << exponent;
			else
				return static_cast<t>(std::pow(T::prime, exponent)); //can be improved!
		}

		//A variant of \ref rotate "rotate" for segments of Eigen vectors
		//The reason we do shiftsegment(V,a,b) as opposed to rotate(V.segment(a,b)) is explained in "In which cases do functions taking a plain Matrix or Array argument fail?" of Eigen documentation. The trick there results in segmentation faults...
		template<typename Derived, typename T>
		void shiftsegment(Eigen::MatrixBase<Derived>& V, int tracker, T rankvalue)
		{
			for (int i = tracker; i < tracker + rankvalue - 1; i++)
				std::swap(V(i), V(i + 1));
		}
	}


	//Transfer the rank to given level.
	template<typename group_t>
	auto transfer(const typename group_t::rank_t& rank, int level) {
		typename group_t::rank_t transferred(rank.size());
		scalar_t<typename group_t::rank_t> cutoff = implementation_details::intexp<group_t>(group_t::power - level);
		for (decltype(transferred.size()) i = 0; i < transferred.size(); i++) //when Eigen gets iterators we can replace this with auto
			transferred[i] = std::min(rank[i], cutoff);
		return transferred;
	}

	//Transfer the differential to given level. We need the ranks both at the original level and the given level, for both domain and range.
	template<typename group_t>
	auto transfer(const typename group_t::diff_t& diff, const typename group_t::rank_t& domain, typename group_t::rank_t& domain_top, const typename group_t::rank_t& range, typename group_t::rank_t& range_top, int level) {
		typedef typename group_t::diff_t diff_t;
		typedef typename diff_t::StorageIndex ind;
		if (diff.size() == 0)
			return diff;

		auto n = implementation_details::intexp<group_t>(group_t::power - level);

		typedef typename std::conditional<SFINAE::is_Sparse<diff_t>::value, sparse_t<diff_t>, diff_t>::type coltype;
		coltype reduceddiff;
		ind track, tracktransfer, sum, j;
		track = tracktransfer = sum = j = 0;
		coltype transfer(summation(range_top), summation(domain_top));

		if constexpr (SFINAE::is_Sparse<diff_t>::value) {
			std::vector<ind> keep(diff.rows(), -1);
			ind counter = 0;
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
			std::vector<ind> keep;
			keep.reserve(range.size() * n);
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
			std::vector<ind> cols;
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
			transfer.prune(0, 0); //remove 0 elements, do not use pruned() !!!!
		return transfer;
	}

	//Transfer the bottom and get all levels.
	template<typename T, typename group_t>
	Levels<T, group_t>::Levels(T& bottom) {
		level.resize(group_t::power + 1);
		level[0] = std::move(bottom);
		for (int i = 1; i < group_t::power + 1; i++)
			level[i] = transfer<group_t>(level[0], i);
	}


	template<typename group_t>
	chains_t<group_t> transfer(const chains_t<group_t>& C, int level) {
		chains_t<group_t> D;
		D.reserve(C.maxindex() + 1);
		for (const auto& rank : C.rank)
			D.rank.push_back(transfer<group_t>(rank, level));
		D.diff.push_back(C.diff.front());
		for (int i = 1; i < C.diff.size(); i++)
			D.diff.push_back(transfer<group_t>(C.diff[i], C.rank[i], D.rank[i], C.rank[i - 1], D.rank[i - 1], level));
		return D;
	}

	//Transfer Junction to the desired level.
	template<typename group_t>
	junction_t<group_t> transfer(const junction_t<group_t>& J, int level) {
		auto rankIn_level = transfer<group_t>(J.rankIn, level);
		auto rank_level = transfer<group_t>(J.rank, level);
		auto rankOut_level = transfer<group_t>(J.rankOut, level);
		auto diffIn_level = transfer<group_t>(J.diffIn, J.rankIn, rankIn_level, J.rank, rank_level, level);
		auto diffOut_level = transfer<group_t>(J.diffOut, J.rank, rank_level, J.rankOut, rankOut_level, level);
		return Junction<typename group_t::rank_t, typename group_t::diff_t>(rank_level, rankOut_level, rankIn_level, diffOut_level, diffIn_level);
	}

	//Transfer generator to level given the ranks at the original level (domain) and the target level (range).
	template<typename group_t, typename Derived>
	Derived transfer(const Eigen::MatrixBase<Derived>& generator, const typename group_t::rank_t& domain, const typename group_t::rank_t& range) {
		Derived transferred(summation(range));
		int trackdom = 0;
		int trackran = 0;
		for (int i = 0; i < domain.size(); i++)
		{
			if (domain(i) == range(i))
				transferred.segment(trackran, range(i)) = static_cast<scalar_t<Derived>>(group_t::prime) * generator.segment(trackdom, range(i));
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

	//Restrict generator to level given the ranks at the original level (domain) and the target level (range).
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

	//Compute the Weyl group action on a generator given the rank of the group it lives in.
	template<typename rank_t, typename Derived>
	Derived action(const Eigen::MatrixBase<Derived>& generator, const rank_t& rank)
	{
		Derived Weyl = generator;
		int tracker = 0;
		for (int i = 0; i < rank.size(); i++) {
			implementation_details::shiftsegment(Weyl, tracker, rank(i));
			tracker += rank(i);
		}
		return Weyl;
	}


	//The inverse of the restriction function on a generator given the ranks at the original level (domain) and the target level (range).
	//Recall that free Mackey functors have injective restrictions, so we only need the generator to be in the image.
	template<typename rank_t, typename T>
	T invRes(const T& generator, const rank_t& domain, const rank_t& range)
	{
		T unrestricted(summation(range));
		int64_t trackdom = 0, trackran = 0;
		for (decltype(range.size()) i = 0; i < range.size(); i++) {
			unrestricted.segment(trackran, range[i]) = generator.segment(trackdom, range[i]);
			trackdom += domain[i];
			trackran += range[i];
		}
		return unrestricted;
	}


	//Writing the transfer of each generator in terms of the generators in the image.
	template<typename group_t>
	dense_t<typename group_t::rank_t> transfer(const Homology<typename group_t::rank_t, typename group_t::diff_t>& low, const Homology<typename group_t::rank_t, typename group_t::diff_t>& high, const typename group_t::rank_t& rank_low, const typename group_t::rank_t& rank_high) {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

		dense_t<rank_t> Tr(high.Groups.number_of_summands(), low.Groups.number_of_summands());
		for (decltype(low.Generators.cols()) i = 0; i < low.Generators.cols(); i++) {
			gen_t<rank_t, diff_t> generator = low.Generators.col(i);
			Tr.col(i) = high.basis(transfer<group_t>(generator, rank_low, rank_high)).transpose();
		}
		return Tr;
	}

	//Writing the restriction of each generator in terms of the generators in the image.
	template<typename group_t>
	dense_t<typename group_t::rank_t> restriction(const Homology<typename group_t::rank_t, typename group_t::diff_t>& high, const Homology<typename group_t::rank_t, typename group_t::diff_t>& low, const typename group_t::rank_t& rank_high, const typename group_t::rank_t& rank_low) {
		dense_t<typename group_t::rank_t> Res(low.Groups.number_of_summands(), high.Groups.number_of_summands());
		for (decltype(high.Generators.cols()) i = 0; i < high.Generators.cols(); i++) {
			gen_t<typename group_t::rank_t, typename group_t::diff_t> generator = high.Generators.col(i);
			Res.col(i) = low.basis(restriction(generator, rank_high, rank_low)).transpose();
		}
		return Res;
	}

	//Writing the Weyl group action on each generator in terms of the other generators.
	template<typename group_t>
	dense_t<typename group_t::rank_t> action(const Homology<typename group_t::rank_t, typename group_t::diff_t>& H, const typename group_t::rank_t& rank) {
		dense_t<typename group_t::rank_t> Weyl(H.Groups.number_of_summands(), H.Groups.number_of_summands());
		for (decltype(H.Generators.cols()) i = 0; i < H.Generators.cols(); i++) {
			gen_t<typename group_t::rank_t, typename group_t::diff_t> generator = H.Generators.col(i);
			Weyl.col(i) = H.basis(action(generator, rank)).transpose();
		}
		return Weyl;
	}
}
