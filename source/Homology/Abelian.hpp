#pragma once
#include "Utility/General.hpp"

///@file
///@brief Contains the class \ref mackey::AbelianGroup

namespace mackey
{
	///	@brief	A finitely generated abelian group
	template<typename T>
	class AbelianGroup {
	public:

		///	@brief	The group represented as an array of numbers 1,2,3,4... where 1 corresponds to Z and n corresponds to Z/n for n>1 (the group then corresponds to their summand)
		///	@note Throughout we assume that group[i] divides group[j] if i<=j
		T group;

		///Default constructor (trivial group)
		AbelianGroup() = default;

		///Given an array of numbers a[0],a[1],...,a[n] constructs the group G(a[0])+G(a[1])+...+G(a[n]) where G(1)=Z and G(n)=Z/n if n>=1
		///We assume that in the provided array, a[i] divides a[j] if i<=j
		AbelianGroup(const T& a);

		///Returns 1 if group is trivial
		bool istrivial() const;

		///Returns 1 if group is cyclic
		bool iscyclic() const;

		///Returns the i-th summand of the group (unmutable). If the i-th summand is Z then this returns 1 and if it's Z/n it returns n.
		auto operator[](int i) const;

		///Returns the i-th summand of the group (mutable reference version). If the i-th summand is Z then this returns 1 and if it's Z/n it returns n.
		auto& operator[](int i);

		///Returns the number of direct summands of the group
		int number_of_summands() const;

		///Compares for isomorphism
		bool operator == (const AbelianGroup<T>& G) const;

		///Compares for non-isomorphism
		bool operator != (const AbelianGroup<T>& G) const;

		///The direct sum of groups. No sorting!
		AbelianGroup<T> operator+(const AbelianGroup<T>& G) const;

		///The order of the group. Returns 0 if it's infinite
		int order() const;

		///The order of an element in the group, provided as an array that expresses it as a linear combination of the generators. Returns 0 if it's infinite
		int order(const T& element) const;

		///Normalize linear map whose range is this group (reduce n in Z/n to 0 etc.).
		template<typename Scalar>
		void normalize(Eigen::Matrix<Scalar, -1, -1>& linmap) const;

		///Normalizes an element in the group to have minimal amount of signs amongst its multiples
		void normalize(T& element) const;

		///Expresses a as a linear combination of elements b[i] in the group. If a is not in their span, returns empty matrix.
		auto span(const T& a, const std::vector<T>& b) const;

		///Returns k if a=k*b in group. Returns 0 if a is not a multiple of b (so make sure a!=0)
		int isMultiple(const T& a, const T& b) const;

		///Returns all automorphisms of the group together with their inverses. Currently only for p-groups and Z+p-groups (otherwise returns empty)
		std::pair <std::vector<dense_t<T>>, std::vector<dense_t<T>>> all_automorphisms() const;

	private:
		///Normalizes an element so that 0<=element[i]<=group[i] if group[i]!=1
		void mod_normalize(T& element) const;

		///Checks if integer matrix induces an isomorphism assuming assuming our group is p-group
		template<typename S>
		bool check_isomorphism_p_group(int p, const S& a) const;

		///Returns all automorphisms of our group assuming it's a p-group
		std::vector<dense_t<T>> aut_p_group(int p) const;

		///Returns the inverses of all provided isos 
		std::vector<dense_t<T>> inverses(const std::vector<dense_t<T>>& isos) const;

		///Returns all primes up to n (not optimized; n is usually very small)
		static std::vector<int> primes(int n);

		///Finds if this abelian group is actually a p-group and returns the p if so (otherwise returns 0)
		int is_p_group() const;

		template<typename S>
		friend std::ostream& operator<<(std::ostream& os, const AbelianGroup<S>& A);
	};
	
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const AbelianGroup<T>& A);
}
#include "impl/Abelian.ipp"
