#pragma once
#include "Chains.hpp"
#include <Eigen/Sparse>

///@file
///@brief Contains the classes \ref mackey::AMT and \ref mackey::EquivariantAMT

namespace mackey {

	///	@brief	 Performs algebraic Morse theory reduction 
	///	@details Reduces given chain complex to a homotopy equivalent one but ideally smaller, using Algebraic Morse Theory
	///	@warning The constructor performs the computation and modifies the chain complex in place to avoid copies. 
	/// If you need the original chain complex please store it before passing it here.
	template<typename diff_t>
	class AMT {
	public:

		///Returns vector of "change of basis" matrices from the original to the reduced basis
		const auto& original_to_reduced() const;

		///Returns vector of "change of basis" matrices from the reduced to the original basis
		const auto& reduced_to_original() const;

		///Returns the average compression ratio as a double in [0,1]. Ideally as close to 0 as possible.
		double reduction_ratio() const;

		///Constructor that reduces
		AMT(std::vector<diff_t>& A, bool compute_original_to_reduced, bool compute_reduced_to_original);

	protected:

		///Constructor that allows only to resize
		AMT(std::vector<diff_t>& A, bool compute_original_to_reduced, bool compute_reduced_to_original, bool onlyresize);

		typedef typename diff_t::StorageIndex ind; ///<The storage type of the differential matrices (eg size_t)
		typedef typename diff_t::Scalar scalar_t; ///<The scalar type of the differential matrices (eg int or Z2)
		typedef Eigen::SparseMatrix<scalar_t, 1, ind> row_major_t; ///<The type of row major of the differential matrices (eg size_t)

		std::vector<diff_t>& diff; ///<The reduced differential
		std::vector<std::map<ind, ind>> morse; ///<The morse matching
		std::vector<std::vector<ind>> critical; ///<The critical basis elements

		///Performs the reduction and sets the "change of basis" f,g.
		void reduce();

		///Sets the morse matchiing and normalizes diff if normalize=1 (this speeds up the AMT algorithm)
		void find_Morse_matching(int k, bool normalize);

		///Normalize diff if not done in find_Morse_matching
		void normalize(int k);

		///Erase the given morse matchings (eg if they are not equivariant).
		void erase_matchings(int k, const std::vector<std::pair<ind, ind>>& toremove);

	private:
		const bool getf;
		const bool getg;

		size_t original_size; ///<Needed for the ratio
		std::vector<diff_t> f; ///<The original to reduced "change of basis"
		std::vector<diff_t> g; ///<The reduced to original "change of basis"
		std::vector<std::map<ind, ind>> morse_inverse; ///<The inverse of the morse matching
		std::vector<std::map<ind, scalar_t>> normalizing_coefficients; ///<The normalizing coefficients used in the computation
		size_t compute_size() const;
		void resize();
		void set_critical();
		void kill_dead_paths();
		auto zig_zag_differential(int k, ind col) const;
		auto zig_zag_f(int k, ind col) const;
		auto zig_zag_g(int k, ind col) const;
		void set_f(int k);
		void set_g(int k);
		void reduce_diff(int k);

	};

	///	@brief	 Performs algebraic Morse theory reduction preserving equivariance
	///	@details Reduces given equivariant chain complex to a homotopy equivalent one but ideally smaller, using equivariant variant of 
	///			 algebraic Morse theory
	///	@warning The constructor performs the computation and modifies the chain complex in place to avoid copies. 
	/// If you need the original chain complex please store it before passing it here.
	template<typename rank_t, typename diff_t>
	class EquivariantAMT : public AMT<diff_t> {

	public:
		///Constructor that performs the computation
		EquivariantAMT(Chains<rank_t, diff_t>& C, bool compute_original_to_reduced, bool compute_reduced_to_original);

	private:

		using typename AMT<diff_t>::ind;
		using typename AMT<diff_t>::scalar_t;

		Chains<rank_t, diff_t>& C;
		std::vector<std::vector<int64_t>> rank_sums;
		auto find_index(int k, ind element) const;
		auto find_orbit(int k, ind element) const;
		auto find_non_equivariant_matching(int k) const;
		void set_rank_sums();
		void compute_ranks();
		void compute();
	};
}
#include "impl/Morse.ipp"
