#pragma once
#include <algorithm>
#include <vector>
#include <array>

///	@file
///	@brief Contains classes for generating permutations, combinations and a factory for such classes

namespace mackey
{

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief			Prototype for coroutine-like iterators that generate elements such as combinations.
	///	@details		Inherit from this class and define a method update() to get a const iterator. \n
	///					You will also need \c begin() and \c end() methods returning such iterators;
	///					\c end() should always be defined by calling the factory \c end() .\n
	///					Example implementations: \c CombinationGenerator and \c PermutationGenerator
	///	@attention		Probably not thread-safe (depends on specialization's \c update() method).
	///	@todo			Implement as coroutine (C++20)
	///	@tparam	spec_t	Used for compile-time polymorphism (CRTP): set it to be the child class.
	///	@tparam	gen_t	The type of the generated element.
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class spec_t, class gen_t>
	class FactoryGenerator
	{
	public:
		const gen_t &operator*() const;						  ///<Returns the generated element
		bool operator!=(const FactoryGenerator &other) const; ///<Inequality of iterators (used to detect if generation has been completed)
		spec_t &operator++();								  ///<Generates next element
		static spec_t end();								  ///<Terminal iterator
	protected:
		bool completed;	 ///<1 if the iterator is end() i.e. if all elements have been generated
		gen_t generated; ///<The currently generated element
	};

	////////////////////////////////////////////////////////////////////////////////////////
	///	@brief		Generates all permutations on a number of letters
	///	@details	Use with a ranged for loop: If \c v is a \c PermutationGenerator object then in
	///				\code for (const auto& i:v) {...} \endcode the variable \c i will range over all permutations.
	///	@warning	Not thread safe!
	///	@todo		Implement via coroutine (C++20)
	///	@tparam T	The data type of our permutations eg \c int
	//////////////////////////////////////////////////////////////////////////////////
	template <class T>
	class PermutationGenerator
	{
		const T n;

	public:
		/// @brief	Computes total number of permutations
		/// @return Factorial \f$n!\f$ where \f$n\f$ is the number of letters
		size_t size() const;

		///	@brief		Constructor sets up the generator
		///	@param n	The total number of letters
		PermutationGenerator(T n);

		///	@brief		Constant iterator that is used in a ranged for loop to generate the permutations.
		///	@warning	Non constant version is illegal
		class ConstIterator : public FactoryGenerator<PermutationGenerator::ConstIterator, std::vector<T>>
		{
			void update();
			friend class PermutationGenerator;													///<Befriending outer class
			friend class FactoryGenerator<PermutationGenerator::ConstIterator, std::vector<T>>; ///<Befriending parent
		};

		/// @brief	Begin iterator
		/// @return An iterator to the first generated element
		ConstIterator begin() const;

		/// @brief	End iterator
		/// @return An iterator to the end of the generator (equality with this indicates that the generator has completed)
		ConstIterator end() const;
	};

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief		Generates all combinations on a number of letters making a number of choices
	///	@details	Use with a ranged for loop: If \c v is a \c CombinationGenerator object then in
	///				\code for (const auto& i:v) {...} \endcode the variable \c i will range over all combinations.
	///	@warning	Not thread safe!
	///	@todo		Implement via coroutine (C++20)
	///	@tparam  T	The data type of our combinations eg \c int
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class T>
	class CombinationGenerator
	{
		const T total, choices;

	public:
		///	@brief	Computes total number of combinations
		///	@return Binomial \f${n}\choose{k}\f$ where \f$n\f$=total and \f$k\f$=choices
		auto size() const;

		///	@brief			Sets up the generator
		///	@param total	The number of letters
		///	@param choices	The number of choices
		CombinationGenerator(T total, T choices);

		///	@brief		Constant iterator that is used in a ranged for loop to generate the combinations.
		///	@warning	Non \c const version is illegal
		class ConstIterator : public FactoryGenerator<CombinationGenerator::ConstIterator, std::vector<T>>
		{
			void update();
			friend class CombinationGenerator;													///<Befriending outer class
			friend class FactoryGenerator<CombinationGenerator::ConstIterator, std::vector<T>>; ///<Befriending parent
			T total, choices;
		};

		/// @brief	Begin iterator
		/// @return An iterator to the first generated element
		ConstIterator begin() const;

		/// @brief	End iterator
		/// @return An iterator to the end of the generator (equality with this indicates that the generator has completed)
		ConstIterator end() const;
	};

	///	@brief		Returns vector of all permutations on \p n letters
	///	@tparam	T	The value type of the permutations
	///	@param	n	The number of letters
	///	@return		Vector of vectors each of which is \f$[a_0,...,a_{n-1}]\f$ with \f$0\le a_i\le n-1\f$ all distinct
	template <class T>
	std::vector<std::vector<T>> all_permutations(T n);

	///	@brief		Returns vector of all combinations on \p n letters choosing \p m many
	///	@tparam	T	The value type of the combinations
	///	@param	n	The number of letters
	///	@param	m	The number of choices
	///	@return		Vector of vectors each of which is \f$[a_0,...,a_{n-1}]\f$ with \f$0\le a_i\le n-1\f$ all distinct
	template <class T>
	std::vector<std::vector<T>> all_combinations(T n, T m);
}
#include "impl/Generators.ipp"