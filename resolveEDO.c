#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

#include "utils.h"
#include "edo.h"
#include "gs.h"

/*
 * resolveEDO.c
 *
 * Lê:
 *  n
 *  a b
 *  ya yb
 *  p q
 *  linhas com r1 r2 r3 r4 (uma ou mais; até EOF)
 *
 * Para cada conjunto r1..r4:
 *  - gera SL (usando genTridiag)
 *  - exibe a matriz aumentada (prnEDOsl)
 *  - resolve via Gauss-Seidel
 *  - imprime solução (linha), número de iterações, norma L2, tempo (ms)
 *
 * Formatos exigidos:
 *  - reais: 15 casas (FORMAT definido em edo.h)
 *  - tempo: 8 casas decimais (%.8e)
 *
 * Observações: fsesetround(FE_DOWNWARD) é usado para forçar arredondamento para baixo.
 */

int main(void)
{
    // Forçar arredondamento para baixo conforme enunciado
    fesetround(FE_DOWNWARD);

    int n;
    scanf("%d", &n); 

    double a, b;
    scanf("%lf %lf", &a, &b);

    double ya, yb;
    scanf("%lf %lf", &ya, &yb);

    double p, q;
    scanf("%lf %lf", &p, &q);

    // ler linhas restantes com 4 coeficientes até EOF
    double r1, r2, r3, r4;
    while ( (scanf("%lf %lf %lf %lf", &r1, &r2, &r3, &r4)) == 4 ) {

        EDo eq;
        eq.n = n;
        eq.a = a; eq.b = b;
        eq.ya = ya; eq.yb = yb;
        eq.p = p; eq.q = q;
        eq.r1 = r1; eq.r2 = r2; eq.r3 = r3; eq.r4 = r4;

        // imprime SL (matriz aumentada)
        prnEDOsl(&eq);

        // gerar estrutura tridiagonal
        Tridiag *sl = genTridiag(&eq);

        // preparar vetor solução
        real_t *x = (real_t *) calloc(n, sizeof(real_t));

        // medir tempo
        rtime_t t0 = timestamp();

        int iters = 0;
        real_t normL2 = 0.0;
        gauss_seidel_tridiag(sl, x, 100, 1e-5, &iters, &normL2);

        rtime_t t_elapsed = timestamp() - t0;

        // imprimir solução (linha única, n elementos)
        for (int i = 0; i < n; ++i) {
            printf(FORMAT, x[i]);
            if (i != n-1) printf(" ");
        }
        printf("\n");

        // número de iterações (apenas o inteiro em linha separada)
        printf("%d\n", iters);

        // norma L2 (15 casas)
        printf(FORMAT "\n", normL2);

        // tempo em milisegundos com 8 casas decimais
        printf("%.8e\n\n", t_elapsed);

        // liberar
        free(sl->D);
        free(sl->Di);
        free(sl->Ds);
        free(sl->B);
        free(sl);
        free(x);
    }

    return 0;
}
