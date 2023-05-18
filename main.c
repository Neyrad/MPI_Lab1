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
    int rank, size;
    int Nt = (int)(t_max / t_step) + 1;
    int Nx = (int)(x_max / x_step) + 1;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double** matrix = (double**)calloc(Nx, sizeof(double*));
    for (int m = 0; m < Nx; ++m)
        matrix[m] = (double*)calloc(Nt, sizeof(double));

    double start = MPI_Wtime();
    fill_matrix_par(matrix, rank, size);
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    printf("Rank: %d -> Time: %0.3f seconds\n", rank, end - start);

    int left_boarder  = (Nx / size) * rank;
    int right_boarder = (Nx / size) * (rank + 1) - 1;

    if (rank == size-1)
        right_boarder = Nx-1;

    if (rank != 0)
        for (int m = left_boarder; m <= right_boarder; ++m)
            MPI_Send(matrix[m], Nt, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);

    if (rank == 0)
        for (int i = 1; i < size; ++i)
        {
            int left_boarder_i  = (Nx / size) * i;
            int right_boarder_i = (Nx / size) * (i + 1) - 1;

            if (i == size-1)
                right_boarder_i = Nx-1;

            for (int m = left_boarder_i; m <= right_boarder_i; ++m)
                MPI_Recv(matrix[m], Nt, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &status);
        }

    if (rank == 0)  
        write_to_out(matrix, Nt, Nx);

    for (int m = 0; m < Nx; ++m)
        free(matrix[m]);
    free(matrix);

    MPI_Finalize();
}
