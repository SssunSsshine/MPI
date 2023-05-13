#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <random>
#include <iostream>
#include <iomanip>


#define dimension 1

using namespace std;

int calñulate_min(int* a, int* b, int len)
{
    int min = a[0] + b[0];
    for (int i = 1; i < len; i++) {

        int tmp = a[i] + b[i];
        min = min < tmp ? min : tmp;
    }
    return min;
}

void fill_array_random(int* arr, int len)
{
    std::random_device rnd_device;
    std::mt19937 mersenne_engine(rnd_device());
    std::uniform_int_distribution<int> dist(1, 50);
    for (int i = 0; i < len; i++) {
        arr[i] = dist(mersenne_engine);
    }
}

void printArr(int* arr, int len)
{
    for (int i = 0; i < len; i++)
    {
        cout << arr[i] << " ";
    }
    cout << "\n";
}

int main(int argc, char* argv[])
{
    int i, j, n, rank, rank_pred, rank_next, current;
    int dims[dimension], periods[dimension];
    MPI_Comm new_comm;
    MPI_Status st;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Îáíóëÿåì ìàññèâ dims è çàïîëíÿåì ìàññèâ periods äëÿ òîïîëîãèè "êîëüöî" 
    for (i = 0; i < dimension; i++)
    {
        dims[i] = 0;
        periods[i] = 1;
    }

    MPI_Dims_create(n, dimension, dims);
    MPI_Cart_create(MPI_COMM_WORLD, dimension, dims, periods, 0, &new_comm);
    MPI_Cart_shift(new_comm, 0, -1, &rank_pred, &rank_next);

    int* A, * B;
    A = (int*)(malloc(sizeof(int) * n));
    B = (int*)(malloc(sizeof(int) * n));


    fill_array_random(A, n);

    fill_array_random(B, n);

    cout << "------------------------------------\n";
    cout << "A " << rank << "\n";
    printArr(A, n);

    cout << "B " << rank << "\n";
    printArr(B, n);
    int maximum = calñulate_min(A, B, n);
    cout << "minimum 0" << " = " << maximum << " ";

    for (j = 1;j < n;j++)
    {
        MPI_Sendrecv_replace(B, n, MPI_INT, rank_next, 2, rank_pred, 2, new_comm, &st);

        cout << "\n\nA\n";
        printArr(A, n);

        cout << "B\n";
        printArr(B, n);

        current = calñulate_min(A, B, n);
        maximum = maximum < current ? current : maximum;

        cout << "minimum " << j << " = " << current << " ";
        cout << "\n";
    }
    cout << "\n";
    cout << "MAXIMUM = " << maximum << "\n";
    MPI_Comm_free(&new_comm);
    free(A);
    free(B);
    MPI_Finalize();
    return 0;
}