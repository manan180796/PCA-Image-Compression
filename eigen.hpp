#ifndef JAMA_EIG_H
#define JAMA_EIG_H


// #include "tnt_array1d.h"
// #include "tnt_DenseMatrix.h"
// #include "tnt_math_utils.h"

#include "DenseMatrix.hpp"

#include <algorithm>
// for min(), max() below

#include <cmath>
// for abs() below

// using namespace TNT;
// using namespace std;


/**
    Computes eigenvalues and eigenvectors of a real (non-complex)
    matrix.
<P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal. That is,
    the diagonal values of D are the eigenvalues, and
    V*V' = I, where I is the identity matrix.  The columns of V
    represent the eigenvectors in the sense that A*V = V*D.

<P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    a + i*b, in 2-by-2 blocks, (a, b; -b, a).  That is, if the complex
    eigenvalues look like
<pre>
          u + iv     .        .          .      .    .
            .      u - iv     .          .      .    .
            .        .      a + ib       .      .    .
            .        .        .        a - ib   .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
</pre>
        then D looks like
<pre>
            u        v        .          .      .    .
           -v        u        .          .      .    .
            .        .        a          b      .    .
            .        .       -b          a      .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
</pre>
    This keeps V a real matrix in both symmetric and non-symmetric
    cases, and A*V = V*D.



    <p>
    The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon the condition number of V.
   <p>
    (Adapted from JAMA, a Java Matrix Library, developed by jointly
    by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
**/

template <class Real>
class Eigenvalue {
    /** Row and column dimension (square matrix).  */
    int n;

    int issymmetric; /* boolean*/

    /** Arrays for internal storage of eigenvalues. */

    DenseMatrix<Real> d; /* real part */
    DenseMatrix<Real> e; /* img part */

    /** Array for internal storage of eigenvectors. */
    DenseMatrix<Real> V;


    // Symmetric Householder reduction to tridiagonal form.

    void tred2() {
        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.

        for (int j = 0; j < n; j++) {
            d(j, 0) = V(n - 1, j);
        }

        // Householder reduction to tridiagonal form.

        for (int i = n - 1; i > 0; i--) {
            // Scale to avoid under/overflow.

            Real scale = 0.0;
            Real h = 0.0;
            for (int k = 0; k < i; k++) {
                scale = scale + std::abs(d(k, 0));
            }
            if (scale == 0.0) {
                e(i, 0) = d(i - 1, 0);
                for (int j = 0; j < i; j++) {
                    d(j, 0) = V(i - 1, j);
                    V(i, j) = 0.0;
                    V(j, i) = 0.0;
                }
            } else {
                // Generate Householder vector.

                for (int k = 0; k < i; k++) {
                    d(k, 0) /= scale;
                    h += d(k, 0) * d(k, 0);
                }
                Real f = d(i - 1, 0);
                Real g = sqrt(h);
                if (f > 0) {
                    g = -g;
                }
                e(i, 0) = scale * g;
                h = h - f * g;
                d(i - 1, 0) = f - g;
                for (int j = 0; j < i; j++) {
                    e(j, 0) = 0.0;
                }

                // Apply similarity transformation to remaining columns.

                for (int j = 0; j < i; j++) {
                    f = d(j, 0);
                    V(j, i) = f;
                    g = e(j, 0) + V(j, j) * f;
                    for (int k = j + 1; k <= i - 1; k++) {
                        g += V(k, j) * d(k, 0);
                        e(k, 0) += V(k, j) * f;
                    }
                    e(j, 0) = g;
                }
                f = 0.0;
                for (int j = 0; j < i; j++) {
                    e(j, 0) /= h;
                    f += e(j, 0) * d(j, 0);
                }
                Real hh = f / (h + h);
                for (int j = 0; j < i; j++) {
                    e(j, 0) -= hh * d(j, 0);
                }
                for (int j = 0; j < i; j++) {
                    f = d(j, 0);
                    g = e(j, 0);
                    for (int k = j; k <= i - 1; k++) {
                        V(k, j) -= (f * e(k, 0) + g * d(k, 0));
                    }
                    d(j, 0) = V(i - 1, j);
                    V(i, j) = 0.0;
                }
            }
            d(i, 0) = h;
        }

        // Accumulate transformations.

        for (int i = 0; i < n - 1; i++) {
            V(n - 1, i) = V(i, i);
            V(i, i) = 1.0;
            Real h = d(i + 1, 0);
            if (h != 0.0) {
                for (int k = 0; k <= i; k++) {
                    d(k, 0) = V(k, i + 1) / h;
                }
                for (int j = 0; j <= i; j++) {
                    Real g = 0.0;
                    for (int k = 0; k <= i; k++) {
                        g += V(k, i + 1) * V(k, j);
                    }
                    for (int k = 0; k <= i; k++) {
                        V(k, j) -= g * d(k, 0);
                    }
                }
            }
            for (int k = 0; k <= i; k++) {
                V(k, i + 1) = 0.0;
            }
        }
        for (int j = 0; j < n; j++) {
            d(j, 0) = V(n - 1, j);
            V(n - 1, j) = 0.0;
        }
        V(n - 1, n - 1) = 1.0;
        e(0, 0) = 0.0;
    }

    // Symmetric tridiagonal QL algorithm.

    void tql2() {
        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.

        for (int i = 1; i < n; i++) {
            e(i - 1, 0) = e(i, 0);
        }
        e(n - 1, 0) = 0.0;

        Real f = 0.0;
        Real tst1 = 0.0;
        Real eps = pow(2.0, -52.0);
        for (int l = 0; l < n; l++) {
            // Find small subdiagonal element

            tst1 = std::max(tst1, std::abs(d(l, 0)) + std::abs(e(l, 0)));
            int m = l;

            // Original while-loop from Java code
            while (m < n) {
                if (std::abs(e(m, 0)) <= eps * tst1) {
                    break;
                }
                m++;
            }


            // If m == l, d(l) is an eigenvalue,
            // otherwise, iterate.

            if (m > l) {
                int iter = 0;
                do {
                    iter = iter + 1;  // (Could check iteration count here.)

                    // Compute implicit shift

                    Real g = d(l, 0);
                    Real p = (d(l + 1, 0) - g) / (2.0 * e(l, 0));
                    Real r = hypot(p, 1.0);
                    if (p < 0) {
                        r = -r;
                    }
                    d(l, 0) = e(l, 0) / (p + r);
                    d(l + 1, 0) = e(l, 0) * (p + r);
                    Real dl1 = d(l + 1, 0);
                    Real h = g - d(l, 0);
                    for (int i = l + 2; i < n; i++) {
                        d(i, 0) -= h;
                    }
                    f = f + h;

                    // Implicit QL transformation.

                    p = d(m, 0);
                    Real c = 1.0;
                    Real c2 = c;
                    Real c3 = c;
                    Real el1 = e(l + 1, 0);
                    Real s = 0.0;
                    Real s2 = 0.0;
                    for (int i = m - 1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e(i, 0);
                        h = c * p;
                        r = hypot(p, e(i, 0));
                        e(i + 1, 0) = s * r;
                        s = e(i, 0) / r;
                        c = p / r;
                        p = c * d(i, 0) - s * g;
                        d(i + 1, 0) = h + s * (c * g + s * d(i, 0));

                        // Accumulate transformation.

                        for (int k = 0; k < n; k++) {
                            h = V(k, i + 1);
                            V(k, i + 1) = s * V(k, i) + c * h;
                            V(k, i) = c * V(k, i) - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * e(l, 0) / dl1;
                    e(l, 0) = s * p;
                    d(l, 0) = c * p;

                    // Check for convergence.

                } while (abs(e(l, 0)) > eps * tst1);
            }
            d(l, 0) = d(l, 0) + f;
            e(l, 0) = 0.0;
        }

        // Sort eigenvalues and corresponding vectors.

        for (int i = 0; i < n - 1; i++) {
            int k = i;
            Real p = d(i, 0);
            for (int j = i + 1; j < n; j++) {
                if (d(j, 0) < p) {
                    k = j;
                    p = d(j, 0);
                }
            }
            if (k != i) {
                d(k, 0) = d(i, 0);
                d(i, 0) = p;
                for (int j = 0; j < n; j++) {
                    p = V(j, i);
                    V(j, i) = V(j, k);
                    V(j, k) = p;
                }
            }
        }
    }


    // Complex scalar division.

    Real cdivr, cdivi;
    void cdiv(Real xr, Real xi, Real yr, Real yi) {
        Real r, d;
        if (abs(yr) > abs(yi)) {
            r = yi / yr;
            d = yr + r * yi;
            cdivr = (xr + r * xi) / d;
            cdivi = (xi - r * xr) / d;
        } else {
            r = yr / yi;
            d = yi + r * yr;
            cdivr = (r * xr + xi) / d;
            cdivi = (r * xi - xr) / d;
        }
    }


  public:
    /** Check for symmetry, then construct the eigenvalue decomposition
    @param A    Square real (non-complex) matrix
    */

    Eigenvalue(const DenseMatrix<Real> &A)
        : n(A.Rows()), V(n, n, 0), d(n, 1, 0), e(n, 1, 0) {
        issymmetric = 1;

        if (issymmetric) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    V(i, j) = A(i, j);
                }
            }

            // Tridiagonalize.
            tred2();

            // Diagonalize.
            tql2();
        }
    }


    /** Return the eigenvector matrix
    @return     V
    */

    DenseMatrix<Real> getV() { return V; }

    /** Return the real parts of the eigenvalues
    @return     real(diag(D))
    */

    DenseMatrix<Real> getRealEigenvalues() { return d; }

    /** Return the imaginary parts of the eigenvalues
    in parameter e_.
    @pararm e_: new matrix with imaginary parts of the eigenvalues.
    */
    DenseMatrix<Real> getImagEigenvalues() { return e; }


    /**
        Computes the block diagonal eigenvalue matrix.
        If the original matrix A is not symmetric, then the eigenvalue
        matrix D is block diagonal with the real eigenvalues in 1-by-1
        blocks and any complex eigenvalues,
        a + i*b, in 2-by-2 blocks, (a, b; -b, a).  That is, if the complex
        eigenvalues look like
    <pre>
              u + iv     .        .          .      .    .
                .      u - iv     .          .      .    .
                .        .      a + ib       .      .    .
                .        .        .        a - ib   .    .
                .        .        .          .      x    .
                .        .        .          .      .    y
    </pre>
            then D looks like
    <pre>
                u        v        .          .      .    .
               -v        u        .          .      .    .
                .        .        a          b      .    .
                .        .       -b          a      .    .
                .        .        .          .      x    .
                .        .        .          .      .    y
    </pre>
        This keeps V a real matrix in both symmetric and non-symmetric
        cases, and A*V = V*D.
        @param D: upon return, the matrix is filled with the block diagonal
        eigenvalue matrix.

    */
    DenseMatrix<Real> getD() {
        DenseMatrix<Real> D = DenseMatrix<Real>(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D(i, j) = 0.0;
            }
            D(i, i) = d(i, 0);
            if (e(i, 0) > 0) {
                D(i, i + 1) = e(i, 0);
            } else if (e(i, 0) < 0) {
                D(i, i - 1) = e(i, 0);
            }
        }
        return D;
    }
};


#endif
// JAMA_EIG_H