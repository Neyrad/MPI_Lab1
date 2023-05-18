#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mpi.h"


void write_to_out(double** matrix, int Nt, int Nx);
void matrix_dump(FILE* out, double** matrix, int Nt, int Nx);
void line_dump(FILE* out, double* line, int Nt);
void fill_matrix_seq(double** matrix, int Nt, int Nx);
void fill_matrix_par(double** matrix, int rank, int size);


void line_dump(FILE* out, double* line, int Nt)
{
    for (int k = 0; k < Nt; ++k)
    {
        fprintf(out, "%.3f", line[k]);
        if (k + 1 != Nt)
            fprintf(out, ",");
    }
}

void matrix_dump(FILE* out, double** matrix, int Nt, int Nx)
{
    for (int m = 0; m < Nx; ++m)
    {
        line_dump(out, matrix[m], Nt);
        fprintf(out, "\n");
    }
}

void write_to_out(double** matrix, int Nt, int Nx)
{
    FILE* f = fopen("output.csv", "w+");
    assert(f);
    matrix_dump(f, matrix, Nt, Nx);
    fclose(f);
}

void fill_matrix_par(double** matrix, int rank, int size)
{
    int Nt = (int)(t_max / t_step) + 1;
    int Nx = (int)(x_max / x_step) + 1;
    int left_boarder  = (Nx / size) * rank;
    int right_boarder = (Nx / size) * (rank + 1) - 1;

    double courant = t_step / x_step;
    double C1 = (courant - 1) / (courant + 1);
    double C2 = 2 * t_step / (courant + 1);

    MPI_Status status;

    if (rank == size-1)
        right_boarder = Nx-1;

    if (rank == 0)
    {
        left_boarder = 1;
        for (int k = 0; k < Nt; ++k)
            matrix[0][k] = ksi(k * t_step);
    }

    for (int k = 0; k < Nt-1; ++k)
    {
        for (int m = left_boarder; m <= right_boarder; ++m)
        {
            if (k == 0)
                matrix[m][k] = fi(m * x_step);

            if ((rank != 0) && (m == left_boarder))
                MPI_Recv(&matrix[m-1][k], 1, MPI_DOUBLE, rank-1, 5, MPI_COMM_WORLD, &status);

            matrix[m][k+1] = matrix[m-1][k] + (matrix[m-1][k+1] - matrix[m][k]) * C1\
                                         + func((k+0.5)*t_step, (m+0.5)*x_step) * C2;
        }

        if (rank != (size-1))
            MPI_Send(&matrix[right_boarder][k], 1, MPI_DOUBLE, rank+1, 5, MPI_COMM_WORLD);
        
    }
}

void fill_matrix_seq(double** matrix, int Nt, int Nx)
{
    double courant = t_step / x_step;
    double C1 = (courant - 1) / (courant + 1);
    double C2 = 2 * t_step / (courant + 1);

    for (int m = 0; m < Nx; ++m)
        matrix[m][0] = fi(m * x_step);

    for (int k = 0; k < Nt; ++k)
        matrix[0][k] = ksi(k * t_step);

    for (int k = 0; k < Nt - 1; ++k)
        for (int m = 1; m < Nx; ++m)
            matrix[m][k+1] = matrix[m-1][k] + (matrix[m-1][k+1] - matrix[m][k]) * C1\
                                         + func((k+0.5)*t_step, (m+0.5)*x_step) * C2;
}