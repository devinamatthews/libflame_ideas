#include <cmath>
#include <algorithm>
#include <type_traits>

/*
 * Minimal matrix structure. Would have tons of bells
 * and whistles IRL.
 */
template <typename T>
struct matrix
{
    T* ptr = nullptr;
    int m = 0;
    int n = 0;
    int rs = 0;
    int cs = 0;
};

/*
 * The usual enums.
 */
typedef enum {UPPER, LOWER} uplo_t;
typedef enum {LEFT, RIGHT} side_t;
typedef enum {UNIT, NONUNIT} diag_t;
typedef enum {NORMAL, TRANSPOSE, CONJUGATE, ADJOINT} conjtrans_t;
typedef enum {TL_TO_BR, BR_TO_TL, /*TR_TO_BL, BL_TO_TR*/} dir_t;

/*
 * A forward declaration is needed so that it can be used in
 * partition_2x2::continue_with.
 */
template <typename T, dir_t Direction, int BS>
struct partition_3x3;

/*
 * Structure representing a 2x2 partitioning of a matrix. The
 * individual partitions are lazily updated by the accessor
 * functions, and the data relevant for looping and point arithmetic
 * is stored as const members.
 *
 * Repartitioning from a 3x3 is done in the continue_with member
 * function, with the movement of the "thick bar" dictated by the
 * Direction template parameter.
 */
template <typename T, dir_t Direction>
struct partition_2x2
{
    T* ptr;
    const int m, rs, cs;
    matrix<T> _ATL, _ATR,
              _ABL, _ABR;
    int m_TL;

    /*
     * The partitions are initialized to the entire matrix A to get
     * copies of rs and cs. The other parameters will be set by the
     * individual accessor functions.
     */
    partition_2x2(matrix<T>& A)
    : ptr(A.ptr), m(A.m), rs(A.rs), cs(A.cs),
      _ATL(A), _ATR(A),
      _ABL(A), _ABR(A),
      m_TL(Direction == TL_TO_BR ? 0 : A.m) {}

    matrix<T>& ATL()
    {
        _ATL.m = _ATL.n = m_TL;
        return _ATL;
    }

    matrix<T>& ABL()
    {
        _ABL.ptr = ptr + m_TL*rs;
        _ABL.m = m - m_TL;
        _ABL.n = m_TL;
        return _ABL;
    }

    /*
     * This function gets either the BL or TR partition based
     * on the triangular/symmetric storage mode. (Note that the
     * "raw" partitions are returned -- a conditional transpose
     * may be necessary to use this function generically).
     */
    template <uplo_t Uplo>
    matrix<T>& ABL_or_TR()
    {
        return (Uplo == LOWER ? ABL() : ATR());
    }

    matrix<T>& ATR()
    {
        _ATR.ptr = ptr + m_TL*cs;
        _ATR.m = m_TL;
        _ATR.n = m - m_TL;
        return _ATR;
    }

    matrix<T>& ABR()
    {
        _ABR.ptr = ptr + m_TL*rs + m_TL*cs;
        _ABR.m = _ABR.n = m - m_TL;
        return _ABR;
    }

    /*
     * Repartition from 3x3 to 2x2. The mapping is determined by
     * the Direction template parameter.
     */
    template <int BS>
    void continue_with(partition_3x3<T, Direction, BS>& A_3x3)
    {
        if (Direction == TL_TO_BR)
        {
            m_TL = A_3x3.m_00 + A_3x3.m_11;
        }
        else
        {
            m_TL = A_3x3.m_00;
        }
    }
};

/*
 * Structure representing a 3x3 partitioning of a matrix. The
 * individual partitions are lazily updated by the accessor
 * functions, and the data relevant for looping and point arithmetic
 * is stored as const members.
 *
 * Repartitioning from a 2x2 is done in the repartition_down member
 * function, with the movement of the "thick bar" dictated by the
 * Direction template parameter.
 *
 * This is the "base class", which is specialized based on whether or
 * not BS = 1. This is the BS != 1 specialization.
 */
template <typename T, dir_t Direction, int BS>
struct partition_3x3_base
{
    T* ptr;
    const int m, rs, cs;
    int m_00, m_11;

    partition_3x3_base(matrix<T>& A)
    : ptr(A.ptr), m(A.m), rs(A.rs), cs(A.cs),
      m_00(Direction == TL_TO_BR ? 0 : A.m), m_11(0) {}

    /*
     * Repartition from 2x2 to 3x3 where the A11 partition is BSxBS
     * (or smaller for the edge case). The partition mapping is determined
     * by the Direction template parameter.
     */
    void repartition_down(partition_2x2<T, Direction>& A_2x2)
    {
        if (Direction == TL_TO_BR)
        {
            m_00 = A_2x2.m_TL;
            m_11 = (BS == 1 ? 1 : std::min(BS, m - m_00));
        }
        else
        {
            m_11 = (BS == 1 ? 1 : std::min(BS, A_2x2.m_TL));
            m_00 = A_2x2.m_TL - m_11;
        }
    }
};

/*
 * This is the BS = 1 specialization of the base class. The size of
 * A11 is a compile time constant 1 as represented by the m_11 member.
 */
template <typename T, dir_t Direction>
struct partition_3x3_base<T, Direction, 1>
{
    T* ptr;
    const int m, rs, cs;
    int m_00;
    constexpr static int m_11 = 1;

    partition_3x3_base(matrix<T>& A)
    : ptr(A.ptr), m(A.m), rs(A.rs), cs(A.cs),
      m_00(Direction == TL_TO_BR ? 0 : A.m) {}

    /*
     * Repartition from 2x2 to 3x3 where the A11 partition is 1x1. The
     * partition mapping is determined by the Direction template parameter.
     */
    void repartition_down(partition_2x2<T, Direction>& A_2x2)
    {
        if (Direction == TL_TO_BR)
        {
            m_00 = A_2x2.m_TL;
        }
        else
        {
            m_00 = A_2x2.m_TL - m_11;
        }
    }
};

/*
 * This is the "real" class that provides the generic stuff.
 */
template <typename T, dir_t Direction, int BS>
struct partition_3x3 : partition_3x3_base<T, Direction, BS>
{
    typedef partition_3x3_base<T, Direction, BS> base;

    matrix<T> _A00, _A01, _A02,
              _A10, _A11, _A12,
              _A20, _A21, _A22;
    using base::ptr;
    using base::m;
    using base::rs;
    using base::cs;
    using base::m_00;
    using base::m_11;

    /*
     * The partitions are initialized to the entire matrix A to get
     * copies of rs and cs. The other parameters will be set by the
     * individual accessor functions.
     */
    partition_3x3(matrix<T>& A)
    : base(A), _A00(A), _A01(A), _A02(A),
               _A10(A), _A11(A), _A12(A),
               _A20(A), _A21(A), _A22(A) {}

    matrix<T>& A00()
    {
        _A00.m = _A00.n = m_00;
        return _A00;
    }

    matrix<T>& A10()
    {
        _A10.ptr = ptr + m_00*rs;
        _A10.m = m_11;
        _A10.n = m_00;
        return _A10;
    }

    /*
     * See notes on partition_2x2::ABL_or_TR.
     */
    template <uplo_t Uplo>
    matrix<T>& A10_or_01()
    {
        return (Uplo == LOWER ? A10() : A01());
    }

    matrix<T>& A20()
    {
        _A20.ptr = ptr + m_00*rs + m_11*rs;
        _A20.m = m - (m_00+m_11);
        _A20.n = m_00;
        return _A20;
    }

    /*
     * See notes on partition_2x2::ABL_or_TR.
     */
    template <uplo_t Uplo>
    matrix<T>& A20_or_02()
    {
        return (Uplo == LOWER ? A20() : A02());
    }

    matrix<T>& A01()
    {
        _A01.ptr = ptr + m_00*cs;
        _A01.m = m_00;
        _A01.n = m_11;
        return _A01;
    }

    /*
     * This version is used when BS = 1, and returns a reference to
     * the single data element alpha_11.
     */
    template <int _BS=BS>
    typename std::enable_if<_BS==1, T&>::type A11()
    {
        return *(ptr + m_00*rs + m_00*cs);
    }

    /*
     * This version is used when BS != 1, and returns a matrix.
     */
    template <int _BS=BS>
    typename std::enable_if<_BS!=1, matrix<T>&>::type A11()
    {
        _A11.ptr = ptr + m_00*rs + m_00*cs;
        _A11.m = _A11.n = m_11;
        return _A11;
    }

    /*
     * See notes on partition_2x2::ABL_or_TR.
     */
    matrix<T>& A21()
    {
        _A21.ptr = ptr + m_00*(rs+cs) + m_11*rs;
        _A21.m = m - (m_00+m_11);
        _A21.n = m_11;
        return _A21;
    }

    /*
     * See notes on partition_2x2::ABL_or_TR.
     */
    template <uplo_t Uplo>
    matrix<T>& A21_or_12()
    {
        return (Uplo == LOWER ? A21() : A12());
    }

    matrix<T>& A02()
    {
        _A02.ptr = ptr + m_00*cs + m_11*cs;
        _A02.m = m_00;
        _A02.n = m - (m_00+m_11);
        return _A02;
    }

    matrix<T>& A12()
    {
        _A12.ptr = ptr + m_00*(rs+cs) + m_11*cs;
        _A12.m = m_11;
        _A12.n = m - (m_00+m_11);
        return _A12;
    }

    matrix<T>& A22()
    {
        _A22.ptr = ptr + m_00*(rs+cs) + m_11*(rs+cs);
        _A22.m = _A22.n = m - (m_00+m_11);
        return _A22;
    }
};

//
// ... sometimes this is really gemv
//
template <typename T, conjtrans_t Op_A, conjtrans_t Op_B, conjtrans_t Op_C=NORMAL>
void gemm(T alpha, const matrix<T>& A, const matrix<T>& B, T beta, matrix<T>& C);

//
// ... sometimes this is really her
//
template <typename T, uplo_t Uplo, conjtrans_t Op_A>
void herk(T alpha, const matrix<T>& A, T beta, matrix<T>& B);

//
// ... dotc since B is a scalar
//
template <typename T, uplo_t Uplo, conjtrans_t Op_A>
void herk(T alpha, const matrix<T>& A, T beta, T& B);

//
// ... sometimes this is really trsv
//
template <typename T, side_t Side, uplo_t Uplo, conjtrans_t Op_A, diag_t Diag>
void trsm(T alpha, const matrix<T>& A, matrix<T>& B);

//
// ... invscal since A is a scalar
//
template <typename T, side_t Side, uplo_t Uplo, conjtrans_t Op_A, diag_t Diag>
void trsm(T alpha, T A, matrix<T>& B);

/*
 * Base case for scalar A. The base case could instead be something
 * that fits in registers (6x6 etc.) and then the trsm/herk/gemm etc.
 * would have the fused vector versions as the base cases.
 */
template <typename T, int Variant, uplo_t Uplo, int BS=1>
void cholesky(T& A)
{
    A = std::sqrt(A);
    // check HPD condition
}

/*
 * 1 for separate code paths for upper and lower storage, 0 for
 * generic code that handles both.
 *
 * For demonstration purposes only; final code would just pick
 * one and go with it. Both compile to the same assembly.
 */
#define SEPARATE_UPPER_AND_LOWER 0

/*
 * This is the meat. Variant, storage type, and block size are all
 * compile-time constants, although the first two don't stricly have
 * to be.
 */
template <typename T, int Variant, uplo_t Uplo, int BS=1>
void cholesky(matrix<T>& A)
{
    /*
     * Save some convenient shorthand for the combined upper/lower case.
     */
    constexpr      side_t   Side = (Uplo == LOWER ?  RIGHT :    LEFT);
    constexpr conjtrans_t HerkOp = (Uplo == LOWER ? NORMAL : ADJOINT);
    constexpr conjtrans_t GemmOp = (Uplo == LOWER ? NORMAL : ADJOINT);

    /*
     * Initialization (worksheet step 4).
     *
     * The 3x3 partitioning is not usually relevant here,
     * but it has to be initialized here to be in the right scope.
     */
    partition_2x2<T, TL_TO_BR    > A_2x2(A);
    partition_3x3<T, TL_TO_BR, BS> A_3x3(A);

    /*
     * Loop guard (worksheet step 3).
     */
    while (A_2x2.m_TL < A.m)
    {
        /*
         * Repartitioning (worksheet step 5a).
         */
        A_3x3.repartition_down(A_2x2);

        /*
         * Update (worksheet step 8).
         */

        #if SEPARATE_UPPER_AND_LOWER

        if (Uplo == LOWER)
        {
            if (Variant == 1)
            {
                /*
                 * We have to save local references to the partitions
                 * we'll need to avoid unecessary recomputation.
                 *
                 * "auto&" is necessary because A11 may be either a
                 * matrix or a scalar.
                 */
                auto& A00 = A_3x3.A00();
                auto& A10 = A_3x3.A10();
                auto& A11 = A_3x3.A11();

                // A10 = A10 * inv( tril( A00 )' )
                trsm<T, RIGHT, LOWER, ADJOINT, NONUNIT>(T(1), A00, A10);

                // A11 = A11 - A10 * A10'
                herk<T, LOWER, NORMAL>(T(-1), A10, T(1), A11);

                // A11 = chol( A11 )
                cholesky<T, 1, LOWER>(A11);
            }
            else if (Variant == 2)
            {
                auto& A10 = A_3x3.A10();
                auto& A11 = A_3x3.A11();
                auto& A20 = A_3x3.A20();
                auto& A21 = A_3x3.A21();

                // A11 = A11 - A10 * A10'
                herk<T, LOWER, NORMAL>(T(-1), A10, T(1), A11);

                // A21 = A21 - A20 * A10'
                gemm<T, NORMAL, ADJOINT>(T(-1), A20, A10, T(1), A21);

                // A11 = chol( A11 )
                cholesky<T, 2, LOWER>(A11);

                // A21 = A21 * inv( tril( A11 )' )
                trsm<T, RIGHT, LOWER, ADJOINT, NONUNIT>(T(1), A11, A21);
            }
            else /* Variant == 3 */
            {
                auto& A11 = A_3x3.A11();
                auto& A21 = A_3x3.A21();
                auto& A22 = A_3x3.A22();

                // A11 = chol( A11 )
                cholesky<T, 3, LOWER>(A11);

                // A21 = A21 * inv( tril( A11 )' )
                trsm<T, RIGHT, LOWER, ADJOINT, NONUNIT>(T(1), A11, A21);

                // A22 = A22 - A21 * A21'
                herk<T, LOWER, NORMAL>(T(-1), A21, T(1), A22);
            }
        }
        else /* Uplo == UPPER */
        {
            if (Variant == 1)
            {
                auto& A00 = A_3x3.A00();
                auto& A01 = A_3x3.A01();
                auto& A11 = A_3x3.A11();

                // A01 = inv( triu( A00 )' ) * A01
                trsm<T, LEFT, UPPER, ADJOINT, NONUNIT>(T(1), A00, A01);

                // A11 = A11 - A01' * A01
                herk<T, UPPER, ADJOINT>(T(-1), A01, T(1), A11);

                // A11 = chol( A11 )
                cholesky<T, 1, UPPER>(A11);
            }
            else if (Variant == 2)
            {
                auto& A01 = A_3x3.A01();
                auto& A11 = A_3x3.A11();
                auto& A02 = A_3x3.A02();
                auto& A12 = A_3x3.A12();

                // A11 = A11 - A01' * A01
                herk<T, UPPER, ADJOINT>(T(-1), A01, T(1), A11);

                // A12 = A12 - A01' * A02
                gemm<T, ADJOINT, NORMAL>(T(-1), A01, A02, T(1), A12);

                // A11 = chol( A11 )
                cholesky<T, 2, UPPER>(A11);

                // A12 = inv( triu( A11 )' ) * A12
                trsm<T, LEFT, UPPER, ADJOINT, NONUNIT>(T(1), A11, A12);
            }
            else /* Variant == 3 */
            {
                auto& A11 = A_3x3.A11();
                auto& A12 = A_3x3.A12();
                auto& A22 = A_3x3.A22();

                // A11 = chol( A11 )
                cholesky<T, 3, UPPER>(A11);

                // A12 = inv( triu( A11 )' ) * A12
                trsm<T, LEFT, UPPER, ADJOINT, NONUNIT>(T(1), A11, A12);

                // A22 = A22 - A12' * A12
                herk<T, UPPER, ADJOINT>(T(-1), A12, T(1), A22);
            }
        }

        #else /* combined upper and lower */

        if (Variant == 1)
        {
            /*
             * The .template notation is unfortunate but unavoidable.
             */
            auto& A00       = A_3x3.A00();
            auto& A10_or_01 = A_3x3.template A10_or_01<Uplo>();
            auto& A11       = A_3x3.A11();

            // lower: A10 = A10 * inv( tril( A00 )' )
            // upper: A01 =       inv( triu( A00 )' ) * A01
            trsm<T, Side, Uplo, ADJOINT, NONUNIT>(T(1), A00, A10_or_01);

            // lower: A11 = A11 - (A10 ) * (A10 )'
            // upper: A11 = A11 - (A01') * (A01')'
            herk<T, Uplo, HerkOp>(T(-1), A10_or_01, T(1), A11);

            // A11 = chol( A11 )
            cholesky<T, 1, Uplo>(A11);
        }
        else if (Variant == 2)
        {
            auto& A10_or_01 = A_3x3.template A10_or_01<Uplo>();
            auto& A11       = A_3x3.A11();
            auto& A20_or_02 = A_3x3.template A20_or_02<Uplo>();
            auto& A21_or_12 = A_3x3.template A21_or_12<Uplo>();

            // lower: A11 = A11 - (A10 ) * (A10 )'
            // upper: A11 = A11 - (A01') * (A01')'
            herk<T, Uplo, HerkOp>(T(-1), A10_or_01, T(1), A11);

            // lower: A21  = A21  - A20 * A10'
            // upper: A12' = A12' - A02 * A01'
            gemm<T, NORMAL, ADJOINT, GemmOp>
                (T(-1), A20_or_02, A10_or_01, T(1), A21_or_12);

            // A11 = chol( A11 )
            cholesky<T, 2, Uplo>(A11);

            // lower: A21 = A21 * inv( tril( A11 )' )
            // upper: A12 =       inv( triu( A11 )' ) * A12
            trsm<T, Side, Uplo, ADJOINT, NONUNIT>(T(1), A11, A21_or_12);
        }
        else /* Variant == 3 */
        {
            auto& A11       = A_3x3.A11();
            auto& A21_or_12 = A_3x3.template A21_or_12<Uplo>();
            auto& A22       = A_3x3.A22();

            // A11 = chol( A11 )
            cholesky<T, 3, Uplo>(A11);

            // lower: A21 = A21 * inv( tril( A11 )' )
            // upper: A12 = inv( triu( A11 )' ) * A12
            trsm<T, Side, Uplo, ADJOINT, NONUNIT>(T(1), A11, A21_or_12);

            // lower: A22 = A22 - A21 * A21'
            // upper: A22 = A22 - A12' * A12
            herk<T, Uplo, HerkOp>(T(-1), A21_or_12, T(1), A22);
        }

        #endif

        /*
         * Repartitioning (worksheet step 5b).
         */
        A_2x2.continue_with(A_3x3);
    }
}

/*
 * Generate some actual code I can disassemble.
 */
template void cholesky<double, 3, LOWER, 1>(matrix<double>&);
