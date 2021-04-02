#pragma once
#include <vector>
#include "Types/Aliases.hpp"
///@file
///@brief Contains the base class for all Smith implementations. 

namespace mackey {

	namespace {
		///Storing a row/column operation
		template<typename T, typename S>
		struct Row_Column_Operation {
			typedef S storage_t;
			//if scalar=0 then this indicates a swap: start<->end
			//otherwise this indicates an operation: end->end-scalar*start
			T scalar;
			S start;
			S end;
			Row_Column_Operation() = default;
			Row_Column_Operation(T scalar, S start, S end) :scalar(scalar), start(start), end(end) {}
		};

		///Storing the pivot and its row/column
		template<typename scalar, typename storage>
		struct Pivot {
			scalar value;
			storage row, col;
			Pivot() = default;
			template<typename T, bool a, bool b>
			Pivot(IteratorNNZ<T, a, b> it) : value(it.value()), row(it.row()), col(it.col()) {}
		};

		//Depending on whether we do full/partial pivoting or we use dense/sparse matrices, the Smith Normal Form may only conditionally need certain members
		template<typename T, bool partial_pivoting = SFINAE::is_finite_cyclic<scalar_t<T>>::value, bool sparse = SFINAE::is_Sparse<T>::value>
		struct smith_conditional_members;

		template<typename T>
		struct smith_conditional_members<T, 0, 0> {
			std::vector<int64_t> Cnorm;
			std::vector<int64_t> Rnorm;
		};

		template<typename T>
		struct smith_conditional_members<T, 0, 1> {
			std::vector<int64_t> Cnorm;
			std::vector<int64_t> Rnorm;
			char current;
			Eigen::SparseMatrix<scalar_t<T>, 1, storage_t<T>> S_row;
			std::vector<Row_Column_Operation<scalar_t<T>, storage_t<T>>> Pops, Qops;
		};

		template<typename T>
		struct smith_conditional_members<T, 1, 0> {
			bool donerow;
			std::vector<int64_t> Cnorm;
		};

		template<typename T>
		struct smith_conditional_members<T, 1, 1> {
			bool donerow;
			std::vector<Row_Column_Operation<scalar_t<T>, storage_t<T>>> Pops, Qops;
		};
	}

	///////////////////////////////////////////////////////////////
	/// The Smith Normal Form and its coefficient matrices

	/// The Smith_Normal_Form of A is a diagonal matrix S with invertible coefficient matrices P,Q,Pi,Qi s.t. P*A*Q=S and Pi*S*Qi=A
	///
	/// S is stored as a diagonal vector. P,Pi are inverses and so are Q,Qi and they are all computed through row-column elimination
	///
	///_S is the type of the original matrix.
	///_R is _S but row major for best performance (required for sparse)
	///_C is _S but column major for best performance (required for sparse)
	/////////////////////////////////////////////////////////////
	template <typename _S, typename _R, typename _C>
	struct Smith_Normal_Form : smith_conditional_members<_S> {
	public:
		row_vector_t<_S> diagonal;		///< The diagonal of the Smith normal form
		_R P;					///< One of the coefficient matrices (S=P*A*Q)
		_R Qi;					///< One of the coefficient matrices (A=Pi*S*Qi)
		_C Q;					///< One of the coefficient matrices (S=P*A*Q)
		_C Pi;					///< One of the coefficient matrices (A=Pi*S*Qi)

		///////////////////////////////////////////////////////////
		///Constructor computes the Smith Normal Form. 
		//
		///The parameters do_P and do_Q can be set to 0 if the P,Pi and Q,Qi matrices are not desired, respectively.
		///The parameter do_sort determines whether we want the diagonal of the Smith Normal Form to be ordered and increasing excluding any units
		///The parameter do_verify determines whether to do an extra check to make sure the Smith Normal Form was computed correctly. 
		///If do_verify=1 then do_P,do_Q are overriden with 1. Thus do_verify=1 is a lot slower than do_verify=0 and should only be used for debugging SNF algos.
		///////////////////////////////////////////////////////////
		Smith_Normal_Form(const _S& S, bool do_P = 1, bool do_Q = 1, bool do_sort = 1, bool do_verify = 0);

	private:

		void verify(const _S&) const;

		typedef _S S_t;
		typedef _R R_t;
		typedef _C C_t;
		typedef typename S_t::StorageIndex ind;

		S_t S;						///< The Smith Normal Form as a diagonal matrix
		const ind M; 				///< Row dimension of original matrix
		const ind N;				///< Column dimension of original matrix
		const bool doP;			///< Whether we want to compute the P and Pi coefficient matrices
		const bool doQ;			///< Whether we want to compute the Q and Qi coefficient matrices

		Pivot<scalar_t<S_t>, ind> piv;

		static constexpr bool partial_pivoting = SFINAE::is_finite_cyclic<typename S_t::Scalar>::value;
		static constexpr bool sparse = SFINAE::is_Sparse<_S>::value;
		static constexpr bool C_exists = !partial_pivoting || !sparse;
		static constexpr bool R_exists = !partial_pivoting;

		typedef std::conditional_t<sparse && !partial_pivoting, typename Eigen::SparseMatrix<scalar_t<S_t>, 1, ind>, int*> S_row_t;

		int64_t metric(ind i, ind j) const;

		void initialize();
		void initialize_norms();
		bool initial_pivoting_per_loop(ind start);
		bool doneRow(ind start);
		bool doneCol(ind start);

		void eliminateRow(ind start);
		void eliminateCol(ind start);
		template<bool row, bool col, char increase_decrease_zero, typename iter>
		void update_norms(iter it);

		void renderPQ();

		void row_pivot(ind start, bool tail = 0);
		void col_pivot(ind start, bool tail = 0);

		void update(int = 0);

		template<bool onlyrow, bool onlycolumn>
		void find_and_set_pivot(ind start);

		void sorter();
		void compute();
	};

}
#include "impl/Smith.ipp"