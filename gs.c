#include <math.h>
#include <stdlib.h>
#include <fenv.h>

#include "gs.h"

/*
 * Gauss-Seidel para sistema tridiagonal representado por Tridiag:
 * Di[i] = sub-diagonal (a_{i,i-1}), D[i] = diagonal (a_{i,i}), Ds[i] = super-diagonal (a_{i,i+1}),
 * B[i] = RHS.
 *
 * Iteração in-place:
 *   x[i] = (B[i] - Di[i]*x[i-1] - Ds[i]*x[i+1_old]) / D[i]
 *
 * Saída:
 *   x[] já preenchido (tamanho n)
 *   iters = número de iterações realizadas
 *   normL2 = norma L2 do resíduo final
 *
 * Retorna 0 em sucesso.
 */
int gauss_seidel_tridiag(Tridiag *sl, real_t *x, int maxit, real_t tol, int *iters, real_t *normL2)
{
    int n = sl->n;
    real_t *Di = sl->Di;
    real_t *D  = sl->D;
    real_t *Ds = sl->Ds;
    real_t *B  = sl->B;

    // inicializa x com zeros
    for (int i = 0; i < n; ++i) x[i] = 0.0;

    // Garantir arredondamento já foi definido pelo chamador (fesetround)
    int k;
    real_t norm = 0.0;

    for (k = 1; k <= maxit; ++k) {
        for (int i = 0; i < n; ++i) {
            real_t s = B[i];
            if (i > 0) s -= Di[i] * x[i-1];
            if (i < n-1) s -= Ds[i] * x[i+1]; // here x[i+1] is old value (Gauss-Seidel uses newest for lower index, old for upper)
            x[i] = s / D[i];
        }

        // calcula resíduo r = A*x - B e norma L2
        norm = 0.0;
        for (int i = 0; i < n; ++i) {
            real_t r = D[i]*x[i];
            if (i > 0) r += Di[i]*x[i-1];
            if (i < n-1) r += Ds[i]*x[i+1];
            r -= B[i];
            norm += r*r;
        }
        norm = sqrt(norm);

        if (norm <= tol) break;
    }

    *iters = (k <= maxit) ? k : maxit;
    *normL2 = norm;

    return 0;
}
