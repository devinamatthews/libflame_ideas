#include <math>
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

enum {UPPER, LOWER};
enum {LEFT, RIGHT};
enum {UNIT, NONUNIT};
enum {NO_TRANSPOSE, TRANSPOSE, CONJUGATE, CONJUGATE_TRANSPOSE};
enum {TL_TO_BR, BR_TO_TL, /*TR_TO_BL, BL_TO_TR*/};

template <typename T, int Direction, int BS>
struct partition_3x3;

template <typename T, int Direction>
struct partition_2x2
{
    matrix<T>& A;
    matrix<T> _ATL, _ATR,
              _ABL, _ABR;
    int m_TL;

    partition_2x2(matrix<T>& A)
    : A(A), _ATL(A), _ATR(A),
            _ABL(A), _ABR(A),
      m_TL(Direction == TL_TO_BR ? 0 : A.m) {}

    matrix<T>& ATL()
    {
        _ATL.m = _ATL.n = m_TL;
        return _ATL;
    }

    matrix<T>& ABL()
    {
        _ABL.ptr = A.ptr + m_TL*A.rs;
        _ABL.m = A.m - m_TL;
        _ABL.n = m_TL;
        return _ABL;
    }

    matrix<T>& ATR()
    {
        _ATR.ptr = A.ptr + m_TL*A.cs;
        _ATR.m = m_TL;
        _ATR.n = A.m - m_TL;
        return _ATR;
    }

    matrix<T>& ABR()
    {
        _ABR.ptr = A.ptr + m_TL*A.rs + m_TL*A.cs;
        _ABR.m = _ABR.n = A.m - m_TL;
        return _ABR;
    }

    template <int BS>
    void repartition_from(partition_3x3<T, Direction, BS>& A_3x3)
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

template <typename T, int Direction, int BS>
struct partition_3x3
{
    matrix<T>& A;
    matrix<T> _A00, _A01, _A02,
              _A10, _A11, _A12,
              _A20, _A21, _A22;
    int m_00, m_11;

    partition_3x3(matrix<T>& A)
    : A(A), _A00(A), _A01(A), _A02(A),
            _A10(A), _A11(A), _A12(A),
            _A20(A), _A21(A), _A22(A),
      m_00(Direction == TL_TO_BR ? 0 : A.m), m_11(0) {}

    matrix<T>& A00()
    {
        _A00.m = _A00.n = m_00;
        return _A00;
    }

    matrix<T>& A10()
    {
        _A10.ptr = A.ptr + m_00*A.rs;
        _A10.m = m_11;
        _A10.n = m_00;
        return _A10;
    }

    matrix<T>& A20()
    {
        _A20.ptr = A.ptr + (m_00+m_11)*A.rs;
        _A20.m = A.m - (m_00+m_11);
        _A20.n = m_00;
        return _A20;
    }

    matrix<T>& A01()
    {
        _A01.ptr = A.ptr + m_00*A.cs;
        _A01.m = m_00;
        _A01.n = m_11;
        return _A01;
    }

    template <int _BS=BS>
    typename std::enable_if<_BS==1, double&>::type A11()
    {
        return *(A.ptr + m_00*A.rs + m_00*A.cs);
    }

    template <int _BS=BS>
    typename std::enable_if<_BS!=1, matrix<T>&>::type A11()
    {
        _A11.ptr = A.ptr + m_00*A.rs + m_00*A.cs;
        _A11.m = _A11.n = m_11;
        return _A11;
    }

    matrix<T>& A21()
    {
        _A21.ptr = A.ptr + (m_00+m_11)*A.rs + m_00*A.cs;
        _A21.m = A.m - (m_00+m_11);
        _A21.n = m_11;
        return _A21;
    }

    matrix<T>& A02()
    {
        _A02.ptr = A.ptr + (m_00+m_11)*A.cs
        _A02.m = m_00;
        _A02.n = A.m - (m_00+m_11);
        return _A02;
    }

    matrix<T>& A12()
    {
        _A12.ptr = A.ptr + m_00*A.rs + (m_00+m_11)*A.cs;
        _A12.m = m_11;
        _A12.n = A.m - (m_00+m_11);
        return _A12;
    }

    matrix<T>& A22()
    {
        _A22.ptr = A.ptr + (m_00+m_11)*A.rs + (m_00+m_11)*A.cs;
        _A22.m = _A22.n = A.m - (m_00+m_11);
        return _A22;
    }

    void repartition_from(partition_2x2<T, Direction>& A_2x2)
    {
        if (Direction == TL_TO_BR)
        {
            m_00 = A_2x2.m_TL;
            m_11 = std::min(BS, A.m - m_00);
        }
        else
        {
            m_11 = std::min(BS, A_2x2.m_TL);
            m_00 = A_2x2.m_TL - m_11;
        }
    }
};

template <typename T, int Op_A, int Op_B>
void gemm(T alpha, const matrix<T>& A, const matrix<T>& B, T beta, matrix<T>& C)
{
    //...
}

template <typename T, int Uplo, int Op_A>
void herk(T alpha, const matrix<T>& A, T beta, matrix<T>& B)
{
    //...
}

template <typename T, int Side, int Uplo, int Op_A, int Diag>
void trsm(T alpha, const matrix<T>& A, matrix<T>& B)
{
    //...
}

template <typename T, int Side, int Uplo, int Op_A, int Diag>
void trsm(T alpha, T A, matrix<T>& B)
{
    //scale by alpha/op(A)...
}

template <typename T, int Variant, int Uplo, int Direction, int BS=1>
void cholesky(T& A)
{
    A = std::sqrt(A);
}

template <typename T, int Variant, int Uplo, int Direction, int BS=1>
void cholesky(matrix<T>& A)
{
    partition_2x2<T, Direction> A_2x2(A);
    partition_3x3<T, Direction, BS> A_3x3(A);

    // TODO: TL -> FF?
    // TODO: check backwards variants
    while (A_2x2.m_TL < A.m)
    {
        A_3x3.repartition_from(A_2x2);

        if (Uplo == LOWER)
        {
            if (Variant == 1)
            {
                auto& A11 = A_3x3.A11();
                auto& A21 = A_3x3.A21();
                auto& A22 = A_3x3.A22();

                cholesky<T, Variant, Uplo, Direction>
                    (A11);
                trsm<T, RIGHT, Uplo, CONJUGATE_TRANSPOSE, NONUNIT>
                    (T(1), A11, A21);
                herk<T, Uplo, NO_TRANSPOSE>
                    (T(-1), A21, T(1), A22);
            }
            else if (Variant == 2)
            {
                auto& A00 = A_3x3.A00();
                auto& A10 = A_3x3.A10();
                auto& A11 = A_3x3.A11();

                trsm<T, RIGHT, Uplo, CONJUGATE_TRANSPOSE, NONUNIT>
                    (T(1), A00, A10);
                herk<T, Uplo, NO_TRANSPOSE>
                    (T(-1), A10, T(1), A11);
                cholesky<T, Variant, Uplo, Direction>
                    (A11);
            }
            else /* Variant == 3 */
            {
                auto& A10 = A_3x3.A10();
                auto& A11 = A_3x3.A11();
                auto& A20 = A_3x3.A20();
                auto& A21 = A_3x3.A21();

                herk<T, Uplo, NO_TRANSPOSE>
                    (T(-1), A10, T(1), A11);
                cholesky<T, Variant, Uplo, Direction>
                    (A11);
                gemm<T, NO_TRANSPOSE, CONJUGATE_TRANSPOSE>
                    (T(-1), A20, A10, T(1), A21);
                trsm<T, RIGHT, Uplo, CONJUGATE_TRANSPOSE, NONUNIT>
                    (T(1), A11, A21);
            }
        }
        else /* Uplo == UPPER */
        {
            if (Variant == 1)
            {
                //TODO
            }
            else if (Variant == 2)
            {
                //TODO
            }
            else /* Variant == 3 */
            {
                //TODO
            }
        }

        A_2x2.repartition_from(A_3x3);
    }
}
