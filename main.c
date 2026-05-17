#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "eliminacaogauss.h"
#include "fatoracaolu.h"
#include "gaussjacobi.h"
#include "gaussseidel.h"

void preencheMatriz(int d, float matriz[d][d]) {
    int i, j;

    for (i = 0; i < d; i++) {
        for (j = 0; j < d; j++) {
            printf("Digite o valor para [%d][%d]: ", i, j);
            scanf("%f", &matriz[i][j]);
        }
    }
}

void preencheVetorB(int d, float vetorb[d]) {
    int i;
    for (i = 0; i < d; i++) {
        printf("Digite o valor para posição %d do vetor b: ", i);
        scanf("%f", &vetorb[i]);
    }
}

void preencheVetorXInicial (int d, float vetorx[d]) {
    int i;
    for (i = 0; i < d; i++) {
        vetorx[i] = 0;
    }
}

void imprimeMatriz(int d, float matriz[d][d]) {
    int i, j;
    for (i = 0; i < d; i++) {
        for (j = 0; j < d; j++) {
            printf("%.2f  ", matriz[i][j]);
        }
        printf("\n");
    }
}

void imprimeVetorB(int d, float vetorb[d]) {
    int i;
    for (i = 0; i < d; i++) {
        printf("%.2f  ", vetorb[i]);
        printf("\n");
    }
}

int main() {
    int dimensaoMtz = 0;

    printf("Informe a dimensao da matriz: ");
    scanf("%d", &dimensaoMtz);
    float matriz[dimensaoMtz][dimensaoMtz];
    float vetorx[dimensaoMtz];
    float vetorb[dimensaoMtz];
    preencheMatriz(dimensaoMtz, matriz);
    preencheVetorB(dimensaoMtz, vetorb);
    preencheVetorXInicial(dimensaoMtz, vetorx);
    imprimeMatriz(dimensaoMtz, matriz);


}