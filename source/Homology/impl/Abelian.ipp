#pragma once
#include "../Abelian.hpp"

///@file
///@brief Contains the methods for finding the automorphisms of finitely generated abelian groups.

namespace mackey
{
	template<typename T>
	AbelianGroup<T>::AbelianGroup(const T& a) :group(a) {}

	template<typename T>
	bool AbelianGroup<T>::istrivial() const { return (group.size() == 0); }

	template<typename T>
	bool AbelianGroup<T>::iscyclic() const { return (group.size() <= 1); }

	template<typename T>
	auto AbelianGroup<T>::operator[](int i) const { return group[i]; }

	template<typename T>
	auto& AbelianGroup<T>::operator[](int i) { return group[i]; }

	template<typename T>
	int AbelianGroup<T>::number_of_summands() const { return group.size(); }

	template<typename T>
	bool AbelianGroup<T>::operator == (const AbelianGroup<T>& G) const {
		if (group.size() != G.group.size())
			return 0;
		return group == G.group;
	}

	template<typename T>
	bool AbelianGroup<T>::operator != (const AbelianGroup<T>& G) const { return !(*this == G); }

	template<typename T>
	AbelianGroup<T> AbelianGroup<T>::operator+(const AbelianGroup<T>& G) const {
		AbelianGroup<T> sum;
		sum.group = T(number_of_summands() + G.number_of_summands());
		for (int i = 0; i < number_of_summands(); i++)
			sum.group[i] = group[i];
		for (int i = number_of_summands(); i < number_of_summands() + G.number_of_summands(); i++)
			sum.group[i] = G.group[i - number_of_summands()];
		return sum;
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& os, const AbelianGroup<T>& A) {
		os << A.group;
		return os;
	}

	template<typename T>
	int AbelianGroup<T>::order() const {
		if (istrivial())
			return 1; //trivial group
		int sum = 0;
		for (int i = 0; i < group.size(); i++) {
			if (group[i] == 1)
				return 0; //infinite order
			else
				sum += group[i];
		}
		return sum;
	}

	template<typename T>
	int AbelianGroup<T>::order(const T& element) const {
		std::vector<int> fractions;
		fractions.reserve(group.size());
		for (int i = 0; i < group.size(); i++) {
			if (group[i] == 1 && element[i] != 0)
				return 0; //infinite order
			else if (element[i] != 0)
				fractions.push_back(group[i] / std::gcd(element[i], group[i]));
		}
		return lcm(fractions);
	}

	template<typename T>
	template<typename Scalar>
	void AbelianGroup<T>::normalize(Eigen::Matrix<Scalar, -1, -1>& linmap) const {
		for (int j = 0; j < linmap.cols(); j++)
			for (int i = 0; i < linmap.rows(); i++)
				if (group[i] != 1)
					linmap(i, j) = (linmap(i, j) + group[i]) % group[i];
	}

	template<typename T>
	void AbelianGroup<T>::mod_normalize(T& element) const {
		for (int i = 0; i < group.size(); i++) {
			if (group[i] != 1)
				element[i] = (element[i] + group[i]) % group[i];
		}
	}

	template<typename T>
	void AbelianGroup<T>::normalize(T& element)const {
		if (element.isZero())
			return;
		mod_normalize(element);
		if (element.isZero())
			return;
		int size = element.size();
		if (size == 1) {
			if (group[0] == 1 && element[0] < 0) {
				element = -element;
				mod_normalize(element);
				return;
			}
			if (group[0] != 1) {
				element[0] = std::gcd(element[0], group[0]);
				return;
			}
		}
		auto n = order(element);
		if (n == 0) {
			for (int i = 0; i < size; i++) {
				if (group[i] == 1 && element[i] < 0) {
					element = -element;
					mod_normalize(element);
					return;
				}
			}
		}
		if (n <= 2)
			return;
		else {
			auto ideal_element = element;
			for (int i = 0; i < size; i++) {
				if (group[i] != 1 && ideal_element[i] != 0)
					ideal_element[i] = std::gcd(ideal_element[i], group[i]);
			}
			int closest_to_ideal = 1;
			int max_hits = -1;
			for (int i = 1; i < n; i++) {
				if (std::gcd(i, n) == 1) {
					int current_hits = 0;
					for (int j = 0; j < size; j++) {
						if ((i * element[j] - ideal_element[j]) % group[j] == 0)
							current_hits++;
					}
					if (max_hits == -1 || max_hits < current_hits) {
						if (current_hits == size) {
							element = ideal_element;
							return;
						}
						max_hits = current_hits;
						closest_to_ideal = i;
					}
				}
			}
			if (closest_to_ideal == 1)
				return;
			element = closest_to_ideal * element;
			mod_normalize(element);
			return;
		}
	}

	template<typename T>
	auto AbelianGroup<T>::span(const T& a, const std::vector<T>& b) const {
		int count = 0;
		for (decltype(group.size()) i = 0; i < group.size(); i++) {
			if (group[i] != 1)
				count++;
		}
		dense_t<T> relation(group.size(), count + b.size());
		relation.setZero();
		int j = 0;
		for (decltype(group.size()) i = 0; i < group.size(); i++) {
			if (group[i] != 1) {
				relation(j, j) = group[i];
				j++;
			}
		}
		for (int i = 0; i < b.size(); i++)
			relation.col(i + count) = b[i];
		Junction<T, dense_t<T>> J;
		J.rank = T::Constant(1, (int)group.size());
		J.diffIn = relation;
		Homology<T, dense_t<T>> H(J, 1);
		gen_t<T, dense_t<T>> cast_a = a.template cast<scalar_t<gen_t<T, dense_t<T>>>>();
		auto v = H.boundary(cast_a);
		if (v.size() == 0)
			return v;
		decltype(v) z(b.size());
		for (int i = 0; i < b.size(); i++)
			z[i] = v[i + count];
		return z;
	}


	template<typename T>
	int AbelianGroup<T>::isMultiple(const T& a, const T& b) const {
		auto u = span(a, { b });
		if (u.size() == 0)
			return 0;
		return (int)u[0];
	}

	template<typename T>
	template<typename S>
	bool AbelianGroup<T>::check_isomorphism_p_group(int p, const S& a) const {
		//first check if it's an endomorphism
		for (int j = 0; j < group.size(); j++) {
			for (int i = 0; i < group.size(); i++) {
				auto quot = group[i] / group[j];
				if (quot != 0 && a(i, j) % quot != 0)
					return 0;
			}
		}
		return (std::gcd(p, abs(static_cast<int>(round(a.template cast<float>().determinant())))) == 1);
	}

	template<typename T>
	std::vector<dense_t<T>> AbelianGroup<T>::aut_p_group(int p) const {
		std::vector<dense_t<T>> isos;
		std::vector<scalar_t<T>> maxelement;
		maxelement.reserve(group.size() * group.size());
		for (int i = 0; i < group.size(); i++) 
			for (int j = 0; j < group.size(); j++) 
				maxelement.push_back(group[j] - 1);
		std::vector<scalar_t<T>> minelement(group.size() * group.size());

		InterpolatingVectorGenerator<std::vector<scalar_t<T>>> allelements(minelement, maxelement);
		isos.reserve(allelements.size());
		for (const auto& i : allelements) {
			dense_t<T> matrix = Eigen::Map<const dense_t<T>>(i.data(), group.size(), group.size());
			if (check_isomorphism_p_group(p, matrix))
				isos.push_back(matrix);
		}
		return isos;
	}

	template<typename T>
	std::vector<dense_t<T>> AbelianGroup<T>::inverses(const std::vector<dense_t<T>>& isos) const {
		dense_t<T> id = dense_t<T>::Identity(isos[0].rows(), isos[0].rows());
		std::vector<dense_t<T>> inverse(isos.size());
		for (size_t i = 0; i < isos.size(); i++)
			if (inverse[i].size() == 0) {
				for (size_t j = i; j < isos.size(); j++) {
					Eigen::Matrix<float, -1, -1> a, b, c;
					Eigen::Matrix<int, -1, -1> d;
					a = isos[i].template cast<float>();
					b = isos[j].template cast<float>();
					c = a * b;
					d = c.template cast<int>();
					normalize(d);
					dense_t<T> product = d.template cast<typename T::Scalar>();
					if (product == id) {
						inverse[i] = isos[j];
						inverse[j] = isos[i];
						break;
					}
				}
		}
		return inverse;
	}

	template<typename T>
	std::vector<int> AbelianGroup<T>::primes(int n) {
		std::vector<int> list;
		list.push_back(2);
		for (int i = 3; i <= n; i++) {
			bool notprime = 0;
			for (const auto& j : list) {
				if (i % j == 0) {
					notprime = 1;
					break;
				}
			}
			if (!notprime)
				list.push_back(i);
		}
		return list;
	}

	template<typename T>
	int AbelianGroup<T>::is_p_group() const {
		int prime = 0;
		if (group[0] == 1)
			return prime;
		auto list = primes(group[0]);
		for (const auto& i : list) {
			if (group[0] % i == 0) {
				prime = i;
				break;
			}
		}
		for (int i = 0; i < group.size(); i++) {
			int current = prime;
			while (current < group[i])
				current *= prime;
			if (current != group[i])
				return 0;
		}
		return prime;
	}

	template<typename T>
	std::pair <std::vector<dense_t<T>>, std::vector<dense_t<T>>> AbelianGroup<T>::all_automorphisms() const {
		std::vector<dense_t<T>> iso, inv;
		auto size = group.size();
		if (size == 0) { //trivial
			auto id = dense_t<T>::Identity(1, 1);
			iso = inv = { id };
		}
		else {
			if (group.minCoeff() > 1) {//finite
				auto p = is_p_group();
				if (p > 0) {
					iso = aut_p_group(p);
					inv = inverses(iso);
				}
			}
			else if (size == 1) { //Z
				dense_t<T> id = dense_t<T>::Identity(1, 1);
				iso = inv = { id,-id };
			}
			else {
				int countZ, locationZ;
				countZ = locationZ = 0;
				for (int i = 0; i < group.size(); i++) {
					if (group[i] == 1) {
						locationZ = i;
						countZ++;
						if (countZ > 1)
							return std::make_pair(iso, inv);
					}
				}
				//finite+Z
				T pgroup(1, size - 1);
				pgroup << group.head(locationZ), group.tail(size - locationZ - 1);
				AbelianGroup pGroup(pgroup);
				auto p = pGroup.is_p_group();
				if (p > 0) {
					auto smalliso = pGroup.aut_p_group(p);
					T onelower = group - T::Constant(size, 1);

					std::vector<scalar_t<T>> maxelement(onelower.data(), onelower.data() + onelower.size()); //the column corresponding to Z
					std::vector<scalar_t<T>> minelement(size);

					InterpolatingVectorGenerator<std::vector<scalar_t<T>>> allcolumns(minelement, maxelement);

					iso.reserve(smalliso.size() * allcolumns.size());
					for (const auto& i : smalliso) {
						for (const auto& j : allcolumns) {
							dense_t<T> matrix(size, size);
							matrix.topLeftCorner(locationZ, locationZ) = i.topLeftCorner(locationZ, locationZ);
							matrix.row(locationZ).setZero();
							matrix.col(locationZ) = Eigen::Map<const col_vector_t<T>>(j.data(), j.size());
							matrix.bottomRightCorner(size - locationZ - 1, size - locationZ - 1) = i.bottomRightCorner(size - locationZ - 1, size - locationZ - 1);
							matrix(locationZ, locationZ) = 1;
							iso.push_back(matrix);
							matrix(locationZ, locationZ) = -1;
							iso.push_back(matrix);
						}
					}
					inv = inverses(iso);
				}
			}
		}
		return std::make_pair(iso, inv);
	}
}
