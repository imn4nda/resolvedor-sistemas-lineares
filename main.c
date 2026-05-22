#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "eliminacaogauss.h"
#include "fatoracaolu.h"
#include "gaussjacobi.h"
#include "gaussseidel.h"
 
void imprimeMatriz(int d, float matriz[d][d]) {
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            printf("%8.4f  ", matriz[i][j]);
        }
        printf("\n");
    }
}
 
void imprimeVetor(int d, float v[d]) {
    for (int i = 0; i < d; i++) {
        printf("  x[%d] = %.6f\n", i, v[i]);
    }
}
 
void copiaMatriz(int d, float src[d][d], float dst[d][d]) {
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            dst[i][j] = src[i][j];
}
 
void copiaVetor(int d, float src[d], float dst[d]) {
    for (int i = 0; i < d; i++)
        dst[i] = src[i];
}
 
int main(int argc, char *argv[]) {
 
    if (argc < 2) {
        fprintf(stderr, "Uso: %s <arquivo_de_entrada>\n", argv[0]);
        return 1;
    }
 
    FILE *fp = fopen(argv[1], "r");
    if (!fp) {
        fprintf(stderr, "Erro ao abrir o arquivo: %s\n", argv[1]);
        return 1;
    }
 
    
    int qtdSistemas, N;
    float precisao;
    if (fscanf(fp, "%d %d %f", &qtdSistemas, &N, &precisao) != 3) {
        fprintf(stderr, "Erro ao ler cabecalho do arquivo.\n");
        fclose(fp);
        return 1;
    }
    printf("Sistemas: %d  |  Dimensao: %d  |  Precisao: %g\n\n",
           qtdSistemas, N, precisao);
 
    
    float (*A)[N]     = malloc(sizeof(float[N][N]));
    float (*Acopia)[N]= malloc(sizeof(float[N][N]));
    float *b          = malloc(N * sizeof(float));
    float *bcopia     = malloc(N * sizeof(float));
    float *x          = malloc(N * sizeof(float));
 
    if (!A || !Acopia || !b || !bcopia || !x) {
        fprintf(stderr, "Erro de alocacao de memoria.\n");
        fclose(fp);
        return 1;
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (fscanf(fp, "%f", &A[i][j]) != 1) {
                fprintf(stderr, "Erro ao ler a matriz A.\n");
                fclose(fp);
                return 1;
            }
 
    printf("Matriz A:\n");
    imprimeMatriz(N, A);
    printf("\n");

    float (*L)[N] = malloc(sizeof(float[N][N]));
    float (*U)[N] = malloc(sizeof(float[N][N]));
    if (!L || !U) {
        fprintf(stderr, "Erro de alocacao para LU.\n");
        fclose(fp);
        return 1;
    }
    copiaMatriz(N, A, Acopia);
    fatoracaoLU(N, Acopia, L, U);
 
    for (int s = 0; s < qtdSistemas; s++) {
 
        for (int i = 0; i < N; i++)
            if (fscanf(fp, "%f", &b[i]) != 1) {
                fprintf(stderr, "Erro ao ler vetor B do sistema %d.\n", s + 1);
                fclose(fp);
                return 1;
            }
 
        printf("========================================\n");
        printf("Sistema %d\n", s + 1);
        printf("Vetor B:\n");
        for (int i = 0; i < N; i++) printf("  b[%d] = %.4f\n", i, b[i]);
        printf("\n");
 
        clock_t inicio, fim;
        double tempoGasto;
 
        
        copiaMatriz(N, A, Acopia);
        copiaVetor(N, b, bcopia);
        inicio = clock();
        eliminacaoGauss(N, Acopia, bcopia, x);
        fim = clock();
        tempoGasto = (double)(fim - inicio) / CLOCKS_PER_SEC;
        printf("--- Eliminacao de Gauss ---\n");
        imprimeVetor(N, x);
        printf("  Tempo: %.6f s\n\n", tempoGasto);
 
        
        copiaVetor(N, b, bcopia);
        inicio = clock();
        resolverLU(N, L, U, bcopia, x);
        fim = clock();
        tempoGasto = (double)(fim - inicio) / CLOCKS_PER_SEC;
        printf("--- Fatoracao LU ---\n");
        imprimeVetor(N, x);
        printf("  Tempo: %.6f s\n\n", tempoGasto);
 
        inicio = clock();
        gaussJacobi(N, A, b, precisao, x);
        fim = clock();
        tempoGasto = (double)(fim - inicio) / CLOCKS_PER_SEC;
        printf("--- Gauss-Jacobi ---\n");
        imprimeVetor(N, x);
        printf("  Tempo: %.6f s\n\n", tempoGasto);
 
        inicio = clock();
        gaussSeidel(N, A, b, precisao, x);
        fim = clock();
        tempoGasto = (double)(fim - inicio) / CLOCKS_PER_SEC;
        printf("--- Gauss-Seidel ---\n");
        imprimeVetor(N, x);
        printf("  Tempo: %.6f s\n\n", tempoGasto);
    }

    free(A); free(Acopia);
    free(b); free(bcopia);
    free(x);
    free(L); free(U);
    fclose(fp);
 
    return 0;
}