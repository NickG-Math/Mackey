#pragma once
#include "Utility/General.hpp" 
#include "Chains.hpp" 
#include "ChangeBasis.hpp"
///@file
///@brief Contains the construction of box products of Chains.


namespace mackey {

	namespace implementation_details{
		//SFINAE aliases

		template<typename T>
		using is_rank = std::enable_if_t<std::is_same<T, row_vector_t<T>>::value, int>;

		template<typename T>
		using is_vector_rank = std::enable_if_t<std::is_same<T, std::vector<row_vector_t<typename T::value_type>>>::value, int>;

		template<typename T>
		using is_arrow = std::enable_if_t<std::is_same<T, arrow_t<T>>::value, int>;

		template<typename T>
		using is_chains = std::enable_if_t<std::is_same<T, chains_t<T>>::value, int>;

		template<typename T>
		using is_junction = std::enable_if_t<std::is_same<T, junction_t<T>>::value, int>;

		template<typename T, bool>
		struct optional_base;

		template<typename T>
		struct optional_base<T, 0> {};

		template<typename T>
		struct optional_base<T, 1> { T _optional; };
	}

	/// @brief				Tensor product of chain complexes and more
	/// @tparam _output_t	The result type of tensor(). Used to resolve "overloads"
	/// @tparam _optional_t The result type of optional(), as long as it's not void. This is the 
	template<typename _output_t, typename _optional_t = void> 
	class Tensor : implementation_details::optional_base<_optional_t, !std::is_same<_optional_t, void>::value> {
	public:
		/// @brief 	Returns the tensor product
		/// @return The tensor product as a non const ref (so you can move it!)
		_output_t& tensor();

		/// @brief 	Returns the optional parameter (detailed rank)
		/// @return The detailed rank as a non const ref (so you can move it!) if _optional_t is not void; raises a static assert otherwise.
		auto& optional();

		/// @brief 	Returns the tensor product
		/// @return The tensor product as a const ref
		const _output_t& tensor() const;

		/// @brief 	Returns the optional parameter (detailed rank)
		/// @return The detailed rank as a non const ref if _optional_t is not void; raises a static assert otherwise.
		const auto& optional() const;

		/// @brief 				Constructor that tensors A,B and optionally at given location only
		/// @tparam _input_t	The type of inputs A,B
		/// @param  A			The first factor of the tensor product
		/// @param  B			The second factor of the tensor product
		/// @param  i			Optionally only perform the tensor at given location. If i==-1 (default) then i is ignored.
		template<typename _input_t>
		Tensor(const _input_t& A, const _input_t& B, int i= -1);

	private:
		static constexpr bool optional_exists = !std::is_same<_optional_t, void>::value;

		_output_t _tensor;

		//The tensor product of the ranks of A,B is the rank of A tensor B
		//Here: _input_t=rank_t, _output_t=rank_t, _optional_t is void, i=-1
		template<typename _input_t, typename implementation_details::is_rank<_input_t> = 0>
		void tensor(const _input_t&, const _input_t&, int = -1);

		//The tensor product of the ranks of A_*,B_* at location i is the rank of (A_* tensor B_*)_i
		//Here: _input_t=std::vector<rank_t>, _output_t=rank_t, _optional_t is std::vector<int64_t> or void, i= location
		template<typename _input_t, typename implementation_details::is_vector_rank<_input_t> = 0, typename S = _output_t, typename implementation_details::is_rank<S> = 0>
		void tensor(const _input_t&, const _input_t&, int = -1);

		//The tensor product of the ranks of A_*,B_* is the rank of A_* tensor B_*
		//Here: _input_t=std::vector<rank_t>, _output_t=std::vector<rank_t>, _optional_t is void, i=-1
		template<typename _input_t, typename T = _output_t, typename S = _output_t, typename implementation_details::is_vector_rank<T> = 0, typename implementation_details::is_vector_rank<S> = 0>
		void tensor(const _input_t&, const _input_t&, int = -1);

		//The tensor product of the Chains A_*,B_* at location i is the arrow (A_* tensor B_*)_i
		//Here: _input_t=Chains<rank_t,diff_t>, _output_t=Arrow<rank_t,diff_t>, _optional_t is std::vector<int64_t> or void, i= location
		template<typename _input_t, typename T = _output_t, typename implementation_details::is_arrow<T> = 0>
		void tensor(const _input_t&, const _input_t&, int = -1);

		//The tensor product of the Chains A_*,B_* up to location i is the Chains (A_* tensor B_*)_* for *<=i
		//Here: _input_t=Chains<rank_t,diff_t>, _output_t=Chains<rank_t,diff_t>, _optional_t is std::vector<std::vector<int64_t>> or void, i= location or -1 (if the whole tensor is desired)
		template<typename _input_t, typename T = _output_t, typename implementation_details::is_chains<T> = 0>
		void tensor(const _input_t&, const _input_t&, int = -1);

		//The tensor product of the Chains A_*,B_* at location i as the Junction (A_* tensor B_*)_i->(A_* tensor B_*)_{i-1}
		//Here: _input_t=Chains<rank_t,diff_t>, _output_t=Junction<rank_t,diff_t>, _optional_t is std::vector<int64_t> or void, i= location
		template<typename _input_t, typename T = _output_t, typename implementation_details::is_junction<T> = 0>
		void tensor(const _input_t&, const _input_t&, int = -1);

	};
}
#include "impl/Box.ipp"