#include <cmath>
#include <algorithm>
#include <type_traits>

template <typename T>
struct matrix
{
    T* ptr = nullptr;
    int m = 0;
    int n = 0;
    int rs = 0;
    int cs = 0;
};

typedef enum {UPPER, LOWER} uplo_t;
typedef enum {LEFT, RIGHT} side_t;
typedef enum {UNIT, NONUNIT} diag_t;
typedef enum {NORMAL, TRANSPOSE, CONJUGATE, ADJOINT} conjtrans_t;
typedef enum {TL_TO_BR, BR_TO_TL, /*TR_TO_BL, BL_TO_TR*/} dir_t;

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
// ... dotc
//
template <typename T, uplo_t Uplo, conjtrans_t Op_A>
void herk(T alpha, const matrix<T>& A, T beta, T& B);

//
// ... sometimes this is really trsv
//
template <typename T, side_t Side, uplo_t Uplo, conjtrans_t Op_A, diag_t Diag>
void trsm(T alpha, const matrix<T>& A, matrix<T>& B);

//
// really invscal
//
template <typename T, side_t Side, uplo_t Uplo, conjtrans_t Op_A, diag_t Diag>
void trsm(T alpha, T A, matrix<T>& B);

template <typename T, int BS=1>
void cholesky(T& A)
{
    A = std::sqrt(A);
    // check HPD condition
}

template <typename T, int BS=1>
void cholesky(matrix<T>& A)
{
    matrix<T> ATL(A), ATR(A),
              ABL(A), ABR(A);
    matrix<T> A00(A), A01(A), A02(A),
              A10(A), A11(A), A12(A),
              A20(A), A21(A), A22(A);

    const int rs_plus_cs = A.rs + A.cs;
    const int m = A.m;

    A21.ptr = A.ptr + BS*A.rs;
    A22.ptr = A.ptr + BS*rs_plus_cs;

    int m_BR = 0;
    for (int m_TL = 0;m_TL < m;m_TL += BS)
    {
        const int m_loc = std::min(BS, m-m_TL);

        m_BR -= m_loc;
        A11.m = A11.n = A21.n = m_loc;
        A21.m = A22.m = A22.n = m_BR;

        // A11 = chol( A11 )
        cholesky<T>(*A11.ptr);

        // A21 = A21 * inv( tril( A11 )' )
        trsm<T, RIGHT, LOWER, ADJOINT, NONUNIT>(T(1), A11, A21);

        // A22 = A22 - A21 * A21'
        herk<T, LOWER, NORMAL>(T(-1), A21, T(1), A22);

        A11.ptr += m_loc*rs_plus_cs;
        A21.ptr += m_loc*rs_plus_cs;
        A22.ptr += m_loc*rs_plus_cs;
    }
}

template void cholesky<double, 1>(matrix<double>&);
