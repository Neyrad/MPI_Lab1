#include <time.h>
#include "func.c"
#include "write.c"

extern double func(double t, double x);
extern double fi(double x);
extern double ksi(double t);

extern double t_max;
extern double x_max;
extern double t_step;
extern double x_step;

int main(int argc, char** argv)
{
    int Nt = (int)(t_max / t_step) + 1;
    int Nx = (int)(x_max / x_step) + 1;

    double** matrix = (double**)calloc(Nx, sizeof(double*));
    for (int m = 0; m < Nx; ++m)
        matrix[m] = (double*)calloc(Nt, sizeof(double));

    clock_t start, end;
    start=clock();
    fill_matrix_seq(matrix, Nt, Nx);
    end=clock();
    printf("Time: %0.3f seconds\n", ((float)(end-start)) / CLOCKS_PER_SEC);
    write_to_out(matrix, Nt, Nx);

    for (int m = 0; m < Nx; ++m)
        free(matrix[m]);
    free(matrix);
}
