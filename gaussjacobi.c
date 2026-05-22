#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaussjacobi.h"
 
void gaussJacobi(int n, float A[n][n], float b[n], float precisao, float x[n]) {
    float *xAnterior = (float *)malloc(n * sizeof(float));
    if (!xAnterior) {
        fprintf(stderr, "Erro de alocacao de memoria em gaussJacobi\n");
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }
 
    int iteracao = 0;
    float erro;
 
    do {
        for (int i = 0; i < n; i++) {
            xAnterior[i] = x[i];
        }

        for (int i = 0; i < n; i++) {
            float soma = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    soma -= A[i][j] * xAnterior[j];
                }
            }
            x[i] = soma / A[i][i];
        }

        float normaDiff = 0.0f, normaX = 0.0f;
        for (int i = 0; i < n; i++) {
            float diff = fabsf(x[i] - xAnterior[i]);
            if (diff > normaDiff) normaDiff = diff;
            if (fabsf(x[i]) > normaX) normaX = fabsf(x[i]);
        }
        /* Evita divisão por zero */
        erro = (normaX > 1e-10f) ? (normaDiff / normaX) : normaDiff;
 
        iteracao++;
    } while (erro > precisao);
 
    printf("  [Jacobi]  iteracoes: %d\n", iteracao);
 
    free(xAnterior);
}