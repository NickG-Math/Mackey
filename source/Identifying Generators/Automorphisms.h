#pragma once
#include "Types/Aliases.h"
#include "Utility/General.h"

///@file
///@brief Contains the methods for finding the automorphisms of finitely generated abelian groups.

namespace Mackey
{

	///Normalize linear map if the range is finitely generated abelian (reduce n in Z/n to 0 etc.).
	template<typename T, typename S>
	void normalize(Eigen::Matrix<T, -1, -1>& linmap, const Eigen::Matrix<S, 1, -1>& range) {
		for (int j = 0; j < linmap.cols(); j++) {
			for (int i = 0; i < linmap.rows(); i++) {
				if (range[i] != 1)
					linmap(i, j) = (linmap(i, j) + range[i]) % range[i];
			}
		}
	}


	///Normalize linear map if the range is finitely generated abelian (reduce n in Z/n to 0 etc.) allowing negative entries
	template<typename Scalar>
	void normalize_minus(Eigen::Matrix<Scalar, -1, -1>& linmap, const Eigen::Matrix<Scalar, 1, -1>& range, bool positive) {
		for (int j = 0; j < linmap.cols(); j++) {
			for (int i = 0; i < linmap.rows(); i++) {
				if (range[i] > 2 && linmap(i, j) == range[i] - 1)
					linmap(i, j) = -1;
				else if (range[i] != 1)
					linmap(i, j) = linmap(i, j) % range[i];
			}
		}
	}




	///Checks if integer matrix induces an isomorphism pGroup->pGroup for finite abelian p-group
	template<typename T, typename S>
	bool check_isomorphism_p_group(int p, const T& a, const S& pGroup) {
		//first check if it's an endomorphism
		for (int j = 0; j < pGroup.size(); j++) {
			for (int i = 0; i < pGroup.size(); i++) {
				auto quot = pGroup[i] / pGroup[j];
				if (quot != 0 && a(i, j) % quot != 0)
					return 0;
			}
		}
		return (std::gcd(p, abs(static_cast<int>(round(a.template cast<float>().determinant())))) == 1);
	}

	///Returns all automorphisms of finite abelian p-group
	template<typename T>
	std::vector<mat_t<T>> aut_p_group(int p, const T& pGroup) {
		std::vector<mat_t<T>> isos;
		std::vector<Scalar_t<T>> maxelement;
		maxelement.reserve(pGroup.size() * pGroup.size());
		for (int i = 0; i < pGroup.size(); i++) {
			for (int j = 0; j < pGroup.size(); j++) {
				maxelement.push_back(pGroup[j] - 1);
			}
		}
		std::vector<Scalar_t<T>> minelement(pGroup.size() * pGroup.size());
		auto allelements = DegreeConstruction(minelement, maxelement);
		isos.reserve(allelements.size());
		for (const auto& i : allelements) {
			mat_t<T> matrix = Eigen::Map<const mat_t<T>>(i.data(), pGroup.size(), pGroup.size());
			if (check_isomorphism_p_group(p, matrix, pGroup))
				isos.push_back(matrix);
		}
		return isos;
	}


	///Returns the inverses of all the automorphisms Group->Group. 
	template<typename T>
	std::vector<mat_t<T>> inverses(const std::vector<mat_t<T>>& isos, const T& Group) {
		mat_t<T> id = mat_t<T>::Identity(isos[0].rows(), isos[0].rows());
		std::vector<mat_t<T>> inverse(isos.size());
		for (int i = 0; i < isos.size(); i++) {
			if (inverse[i].size() == 0) {
				for (int j = i; j < isos.size(); j++) {
					Eigen::Matrix<float, -1, -1> a, b, c;
					Eigen::Matrix<int, -1, -1> d;
					a = isos[i].template cast<float>();
					b = isos[j].template cast<float>();
					c = a * b;
					d = c.template cast<int>();
					normalize(d, Group);
					mat_t<T> product = d.template cast<typename T::Scalar>();
					if (product == id) {
						inverse[i] = isos[j];
						inverse[j] = isos[i];
						break;
					}
				}
			}
		}
		return inverse;
	}



	///Returns all primes up to n (not optimized; n is usually very small)
	std::vector<int> primes(int n) {
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

	///Finds if given abelian group is actually a p-group and returns the p if so (otherwise returns 0)
	template<typename T>
	int is_p_group(const T& Group) {
		int prime = 0;
		if (Group[0] == 1)
			return prime;
		auto list = primes(Group[0]);
		for (const auto& i : list) {
			if (Group[0] % i == 0) {
				prime = i;
				break;
			}
		}
		for (int i = 0; i < Group.size(); i++) {
			int current = prime;
			while (current < Group[i]) {
				current *= prime;
			}
			if (current != Group[i])
				return 0;
		}
		return prime;
	}





	///Returns all automorphisms of a finitely generated abelian group together with their inverses. 
	///Currently only for p-groups and Z+p-groups (otherwise returns empty)
	template<typename T>
	std::pair <std::vector<mat_t<T>>, std::vector<mat_t<T>>> all_automorphisms(const T& Group) {
		std::vector<mat_t<T>> iso, inv;
		auto size = Group.size();
		if (size == 0) { //trivial
			auto id = mat_t<T>::Identity(1, 1);
			iso = inv = { id };
		}
		else {
			if (Group.minCoeff() > 1) {//finite
				auto p = is_p_group(Group);
				if (p > 0) {
					iso = aut_p_group(p, Group);
					inv = inverses(iso, Group);
				}
			}
			else if (size == 1) { //Z
				mat_t<T> id = mat_t<T>::Identity(1, 1);
				iso = inv = { id,-id };
			}
			else {
				int countZ, locationZ;
				countZ = locationZ = 0;
				for (int i = 0; i < Group.size(); i++) {
					if (Group[i] == 1) {
						locationZ = i;
						countZ++;
						if (countZ > 1)
							return std::make_pair(iso, inv);
					}
				}
				//finite+Z
				T pGroup(1, size - 1);
				pGroup << Group.head(locationZ), Group.tail(size - locationZ - 1);
				auto p = is_p_group(pGroup);
				if (p > 0) {
					auto smalliso = aut_p_group(p, pGroup);
					T onelower = Group - T::Constant(size, 1);

					std::vector<Scalar_t<T>> maxelement(onelower.data(), onelower.data() + onelower.size()); //the column corresponding to Z
					std::vector<Scalar_t<T>> minelement(size);
					auto allcolumns = DegreeConstruction(minelement, maxelement);

					iso.reserve(smalliso.size() * allcolumns.size());
					for (const auto& i : smalliso) {
						for (const auto& j : allcolumns) {
							mat_t<T> matrix(size, size);
							matrix.topLeftCorner(locationZ, locationZ) = i.topLeftCorner(locationZ, locationZ);
							matrix.row(locationZ).setZero();
							matrix.col(locationZ) = Eigen::Map<const col_t<T>>(j.data(), j.size());
							matrix.bottomRightCorner(size - locationZ - 1, size - locationZ - 1) = i.bottomRightCorner(size - locationZ - 1, size - locationZ - 1);
							matrix(locationZ, locationZ) = 1;
							iso.push_back(matrix);
							matrix(locationZ, locationZ) = -1;
							iso.push_back(matrix);
						}
					}
					inv = inverses(iso, Group);
				}
			}
		}
		return std::make_pair(iso, inv);
	}
}
