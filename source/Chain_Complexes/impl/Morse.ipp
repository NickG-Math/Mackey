#pragma once
#include "../Morse.hpp"

namespace mackey {

	namespace implementation_details {
		template<typename T>
		inline bool unit(T a) {
			if constexpr (std::is_integral_v<T>)
				return (a == 1 || a == -1);
			else if constexpr (SFINAE::is_finite_cyclic<T>::value || std::is_floating_point_v<T>)
				return (a != 0);
			static_assert(std::is_integral_v<T> || SFINAE::is_finite_cyclic<T>::value, "Currently only supported for Z and field coefficients (Z/N, floats)");
		}

		template<typename T, typename ind>
		ind location(const std::vector<T>& object, int k, ind i) {
			if (k < 0 || k >= object.size())
				return -1;
			if constexpr (std::is_same_v<T, std::vector<ind>>) {
				auto it = std::lower_bound(object[k].begin(), object[k].end(), i);
				if (it != object[k].end() && *it == i)
					return (it - object[k].begin());
				return ind(-1);
			}
			else {
				auto it = object[k].find(i);
				if (it != object[k].end())
					return it->second;
				return -1;
			}
		}
	}


	template<typename diff_t>
	const auto& AMT< diff_t>::original_to_reduced() const {
		if (!getf) {
			std::cerr << "You said you didn't want the original to reduced!";
			abort();
		}
		return f;
	}

	template<typename diff_t>
	const auto& AMT< diff_t>::reduced_to_original() const {
		if (!getg) {
			std::cerr << "You said you didn't want the reduced to original!";
			abort();
		}
		return g;
	}

	template<typename diff_t>
	double AMT< diff_t>::reduction_ratio() const {
		auto reduced_size = compute_size();
		return (double)reduced_size / (double)original_size;
	}

	template<typename diff_t>
	AMT< diff_t>::AMT(std::vector<diff_t>& A, bool compute_original_to_reduced, bool compute_reduced_to_original, bool onlyresize) : diff(A), getf(compute_original_to_reduced), getg(compute_reduced_to_original) {
		if (diff.size() == 1)
			return;
		resize();
		if (onlyresize)
			return;
		for (int i = 1; i < diff.size(); i++)
			find_Morse_matching(i, 1); //find the matching and normalize
		reduce(); //do the reduction using what you found
	}

	template<typename diff_t>
	AMT< diff_t>::AMT(std::vector<diff_t>& A, bool compute_original_to_reduced, bool compute_reduced_to_original) : AMT(A, compute_original_to_reduced, compute_reduced_to_original, 0) {}

	template<typename diff_t>
	void AMT< diff_t>::reduce() {
		set_critical();
		kill_dead_paths();
		if (getf) {
			f.resize(diff.size());
			for (int k = 0; k < diff.size(); k++)
				set_f(k);
		}
		if (getg) {
			g.resize(diff.size());
			for (int k = 0; k < diff.size(); k++)
				set_g(k);
		}
		for (int k = 1; k < diff.size(); k++)
			reduce_diff(k);
	}

	template<typename diff_t>
	void AMT< diff_t>::find_Morse_matching(int k, bool normalize) {
		static_assert(SFINAE::is_Sparse<diff_t>::value, "AMT currently only supported for sparse matrices!");
		row_major_t diff_row = diff[k];
		for (ind col = 0; col < diff[k].cols(); col++) { //for each column
			if (diff[k].col(col).nonZeros() == 0)
				continue;
			typename diff_t::ReverseInnerIterator it(diff[k], col); //last in column
			typename row_major_t::InnerIterator it_r(diff_row, it.row()); //first in row
			if (it_r.col() == col && implementation_details::unit(it_r.value())) { //first in its row and last in its column, also unit
				morse[k - 1][it.col()] = it.row();
				morse_inverse[k - 1][it.row()] = it.col();
				if (it.value() != -1) {//normalize the weight
					diff_row.row(it.row()) *= (scalar_t(-1) / it.value());
					normalizing_coefficients[k - 1][it.row()] = scalar_t(-1) / it.value();
				}
				diff_row.coeffRef(it.row(), it.col()) = 0; //remove entry
			}
		}
		diff_row.prune(0, 0);
		if (normalize)
			diff[k] = diff_row;
	}

	template<typename diff_t>
	void AMT< diff_t>::normalize(int k) {
		row_major_t diff_row = diff[k];
		for (auto it = morse[k - 1].begin(); it != morse[k - 1].end(); it++) {
			auto itnorm = normalizing_coefficients[k - 1].find(it->second); //normalizing coefficients ONLY stores the coefficients that are NOT 1.
			if (itnorm != normalizing_coefficients[k - 1].end())
				diff_row.row(it->second) *= itnorm->second;
			diff_row.coeffRef(it->second, it->first) = 0;
		}
		diff_row.prune(0, 0);
		diff[k] = diff_row;
	}

	template<typename diff_t>
	void AMT< diff_t>::erase_matchings(int k, const std::vector<std::pair<ind, ind>>& toremove) {
		for (auto i : toremove) {
			morse[k - 1].erase(i.first);
			morse_inverse[k - 1].erase(i.second);
			normalizing_coefficients[k - 1].erase(i.second);
		}
	}

	template<typename diff_t>
	size_t AMT< diff_t>::compute_size() const {
		size_t size = 0;
		size += diff[1].rows();
		for (const auto& i : diff)
			size += i.cols();
		return size;
	}

	template<typename diff_t>
	void AMT< diff_t>::resize() {
		morse.resize(diff.size() - 1);
		morse_inverse.resize(diff.size() - 1);
		normalizing_coefficients.resize(diff.size() - 1);
		critical.resize(diff.size());
		original_size = compute_size();
	}



	template<typename diff_t>
	void AMT< diff_t>::set_critical() {
		for (int i = 0; i < diff.size(); i++) {
			auto size = (i == 0) ? diff[1].rows() : diff[i].cols();
			critical[i].reserve(size);
			for (ind j = 0; j < size; j++) {
				if (implementation_details::location(morse_inverse, i, j) == -1 && implementation_details::location(morse, i - 1, j) == -1)
					critical[i].push_back(j);
			}
		}
	}

	template<typename diff_t>
	void AMT< diff_t>::kill_dead_paths() {
		//set rows/columns to 0 if they lead to dead end paths
		for (int i = 1; i < diff.size(); i++)
			diff[i].prune([&](ind row, ind col, scalar_t)
				{
					return (implementation_details::location(morse, i - 2, row) == -1 && implementation_details::location(morse_inverse, i, col) == -1);
				}
		);
	}

	template<typename diff_t>
	auto AMT< diff_t>::zig_zag_differential(int k, ind col) const { //computes all paths starting with fixed v in I_k' and ending in I_{k-1}'
		std::map<ind, scalar_t> final_paths; //we encode a path from v to u of total weight w through the pair (u,w)
		std::vector<std::pair<ind, scalar_t>> temp_paths; //we encode a path from v to u of total weight w as paths[v][u]=w
		temp_paths.reserve(diff[k].rows());
		temp_paths.emplace_back(col, 1); //initial path
		while (!temp_paths.empty()) { //find everything it's connected to
			auto pair = temp_paths.back();
			temp_paths.pop_back();
			for (typename diff_t::InnerIterator it(diff[k], pair.first); it; ++it) {
				auto next_pair = std::make_pair(it.row(), it.value());
				auto loc = implementation_details::location(morse_inverse, k - 1, next_pair.first);
				if (loc == -1)  //path terminates because critical
					final_paths[next_pair.first] += pair.second * next_pair.second;
				else  //not critical, go one step further in this path
					temp_paths.emplace_back(loc, pair.second * next_pair.second);
			}
		}
		return final_paths;
	}

	template<typename diff_t>
	auto AMT< diff_t>::zig_zag_f(int k, ind col) const { //computes all paths starting with fixed v in I_k' and ending in I_k
		std::map<ind, scalar_t> final_paths; //we encode a path from v to u of total weight w through the pair (u,w)
		std::vector<std::pair<ind, scalar_t>> temp_paths; //we encode a path from v to u of total weight w as paths[v][u]=w
		temp_paths.reserve(diff[k].rows());
		temp_paths.emplace_back(col, 1); //initial path
		final_paths[col] = 1; //initial path
		while (!temp_paths.empty()) { //find everything it's connected to
			auto pair = temp_paths.back();
			temp_paths.pop_back();
			for (typename diff_t::InnerIterator it(diff[k], pair.first); it; ++it) {
				auto next_pair = std::make_pair(it.row(), it.value());
				auto loc = implementation_details::location(morse_inverse, k - 1, next_pair.first);
				if (loc != -1) { //not critical, go one step further in this path
					temp_paths.emplace_back(loc, pair.second * next_pair.second);
					final_paths[loc] += pair.second * next_pair.second;
				}
			}
		}
		return final_paths;
	}


	template<typename diff_t>
	auto AMT< diff_t>::zig_zag_g(int k, ind col) const { //computes all paths starting with fixed v in I_k and ending in I_k'
		std::map<ind, scalar_t> final_paths; //we encode a path from v to u of total weight w through the pair (u,w)
		std::vector<std::pair<ind, scalar_t>> temp_paths; //we encode a path from v to u of total weight w as paths[v][u]=w
		temp_paths.reserve(diff[k].rows());
		temp_paths.emplace_back(col, 1); //initial path

		while (!temp_paths.empty()) { //exploration
			auto pair = temp_paths.back();
			temp_paths.pop_back();
			//first is inverse match
			auto loc = implementation_details::location(morse_inverse, k, pair.first);
			if (loc != -1) { //go one step back
				for (typename diff_t::InnerIterator it(diff[k + 1], loc); it; ++it) {
					auto next_pair = std::make_pair(it.row(), it.value());
					temp_paths.emplace_back(next_pair.first, pair.second * next_pair.second);
				}
			}
			else { //check if critical
				auto loccrit = implementation_details::location(critical, k, pair.first);
				if (loccrit != -1)
					final_paths[loccrit] += pair.second;
			}
		}
		return final_paths;
	}

	template<typename diff_t>
	void AMT< diff_t>::set_f(int k) {
		std::vector<Eigen::Triplet<scalar_t, ind>> triplets;
		triplets.reserve(2 * critical[k].size());
		if (k == 0) {
			f[k].resize(diff[1].rows(), critical[0].size());
			for (ind index = 0; index < critical[0].size(); index++)
				triplets.emplace_back(critical[0][index], index, 1);
		}
		else {
			f[k].resize(diff[k].cols(), critical[k].size());
			for (ind index = 0; index < critical[k].size(); index++) {
				auto paths = zig_zag_f(k, critical[k][index]);
				for (const auto& pairs : paths) {
					if (pairs.second != 0)
						triplets.emplace_back(pairs.first, index, pairs.second);
				}
			}
		}
		f[k].setFromTriplets(triplets.begin(), triplets.end());
	}


	template<typename diff_t>
	void AMT< diff_t>::set_g(int k) {
		std::vector<Eigen::Triplet<scalar_t, ind>> triplets;
		triplets.reserve(2 * critical[k].size());
		if (k < diff.size() - 1) {
			g[k].resize(critical[k].size(), diff[k + 1].rows());
			for (ind j = 0; j < diff[k + 1].rows(); j++) {
				auto paths = zig_zag_g(k, j);
				for (const auto& pairs : paths) {
					if (pairs.second != 0) {
						if (normalizing_coefficients[k].find(j) == normalizing_coefficients[k].end())
							triplets.emplace_back(pairs.first, j, pairs.second);
						else
							triplets.emplace_back(pairs.first, j, normalizing_coefficients[k][j] * pairs.second);
					}
				}
			}
		}
		else {
			g.back().resize(critical.back().size(), diff.back().cols());
			for (ind index = 0; index < critical.back().size(); index++)
				triplets.emplace_back(index, critical.back()[index], 1);
		}
		g[k].setFromTriplets(triplets.begin(), triplets.end());
	}


	template<typename diff_t>
	void AMT< diff_t>::reduce_diff(int k) {
		std::vector<Eigen::Triplet<scalar_t, ind>> triplets;
		triplets.reserve(diff[k].nonZeros());
		for (ind index = 0; index < critical[k].size(); index++) {
			auto paths = zig_zag_differential(k, critical[k][index]);
			for (const auto& pairs : paths) {
				if (pairs.second != 0)
					triplets.emplace_back(implementation_details::location(critical, k - 1, pairs.first), index, pairs.second);
			}
		}
		diff[k].resize(critical[k - 1].size(), critical[k].size());
		diff[k].setFromTriplets(triplets.begin(), triplets.end());
	}

	template<typename rank_t, typename diff_t>
	EquivariantAMT< rank_t, diff_t>::EquivariantAMT(Chains<rank_t, diff_t>& C, bool compute_original_to_reduced, bool compute_reduced_to_original)
		: AMT<diff_t>(C.diff, compute_original_to_reduced, compute_reduced_to_original, 1), C(C) {
		if (this->diff.size() == 1)
			return;
		compute();
	}

	template<typename rank_t, typename diff_t>
	auto EquivariantAMT< rank_t, diff_t>::find_index(int k, ind element) const {
		return std::upper_bound(rank_sums[k].begin(), rank_sums[k].end(), element) - 1 - rank_sums[k].begin();
	}

	template<typename rank_t, typename diff_t>
	auto EquivariantAMT< rank_t, diff_t>::find_orbit(int k, ind element) const {
		auto small_index = find_index(k, element);
		std::vector<ind> orbit;
		auto size = C.rank[k][small_index];
		orbit.reserve(size);
		for (size_t i = 0; i < size; i++) {
			orbit.push_back(element);
			element++;
			if (element >= rank_sums[k][small_index + 1])
				element = rank_sums[k][small_index];
		}
		return orbit;
	}

	template<typename rank_t, typename diff_t>
	auto EquivariantAMT< rank_t, diff_t>::find_non_equivariant_matching(int k) const {
		std::vector<std::pair<ind, ind>> must_remove_key_value_pairs;
		must_remove_key_value_pairs.reserve(this->morse[k - 1].size());
		for (auto it = this->morse[k - 1].begin(); it != this->morse[k - 1].end(); it++) { //go through to see which are equivariant and which aren't
			auto orbit_domain = find_orbit(k, it->first);
			auto orbit_range = find_orbit(k - 1, it->second);
			if (orbit_domain.size() != orbit_range.size())
				must_remove_key_value_pairs.emplace_back(it->first, it->second);
			else
				for (size_t i = 0; i < orbit_domain.size(); i++) {
					auto it_find = this->morse[k - 1].find(orbit_domain[i]);
					if (it_find == this->morse[k - 1].end() || it_find->second != orbit_range[i]) {
						must_remove_key_value_pairs.emplace_back(it->first, it->second);
						break;
					}
				}
		}
		return must_remove_key_value_pairs;
	}

	template<typename rank_t, typename diff_t>
	void EquivariantAMT< rank_t, diff_t>::set_rank_sums() {
		rank_sums.reserve(C.rank.size());
		for (const auto& i : C.rank) {
			std::vector<int64_t> curr(i.size() + 1);
			int64_t sum = 0;
			for (int s = 0; s < i.size(); s++) {
				curr[s] = sum;
				sum += i[s];
			}
			curr[i.size()] = sum;
			rank_sums.push_back(curr);
		}
	}


	template<typename rank_t, typename diff_t>
	void EquivariantAMT< rank_t, diff_t>::compute_ranks() {
		for (int k = 0; k < C.rank.size(); k++) {
			rank_t rankhere(this->critical[k].size());
			size_t counter = 0;
			for (size_t i = 0; i < this->critical[k].size();) {
				auto index = find_index(k, this->critical[k][i]);
				rankhere[counter] = C.rank[k][index];
				i += C.rank[k][index];
				counter++;
			}
			rankhere.conservativeResize(counter);
			C.rank[k] = rankhere;
		}
	}

	template<typename rank_t, typename diff_t>
	void EquivariantAMT< rank_t, diff_t>::compute() {
		set_rank_sums(); //initial setup
		for (int k = 1; k < C.diff.size(); k++) {
			this->find_Morse_matching(k, 0); //nonequivariant matching, dont normalize
			this->erase_matchings(k, find_non_equivariant_matching(k));; //find the non equivariant matchings and erase them
			this->normalize(k); //normalize matrices
		}
		this->reduce(); //reduce
		compute_ranks();	//now compute the ranks
	}
}