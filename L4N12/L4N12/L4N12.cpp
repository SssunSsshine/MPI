#include "mpi.h"
#include <iostream>
#include <random>
#include <iomanip>

using namespace std;

int main(int argc, char** argv)
{
    const int n = 5;
    if (n <= 2) {
        return 0;
    }
    const int m = n + n - 1;
    float a[n][n], b[m];
    int xpose, rank;
    int array_of_displacements[n], array_of_blocklengths[n];

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = j + 10 * i;
        }
    }

    if (rank == 0) {
        cout << "A\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << left
                    << setw(5)
                    << a[i][j]
                    << " ";
            }
            cout << "\n\n";
        }
    }

    int mid = n / 2;
    int add_num = 1;
    int num_index = n + 1;
    if (n % 2 == 0) {
        --mid;
        add_num = 2;
        num_index = n;
  
    }

    array_of_blocklengths[0] = 1;
    array_of_displacements[0] = mid;
    for (int i = 1; i < n; i++) {
        array_of_blocklengths[i] = 1;
        array_of_displacements[i] = array_of_displacements[i - 1] + 2 * mid + add_num;
    }
    array_of_blocklengths[num_index/ 2 - 1] = n;
    array_of_displacements[num_index / 2 - 1] = array_of_displacements[num_index / 2 - 2] + mid + add_num;

    MPI_Type_indexed(n, array_of_blocklengths, array_of_displacements, MPI_FLOAT, &xpose);
    MPI_Type_commit(&xpose);
    MPI_Sendrecv(a, 1, xpose, rank, 0, b, m, MPI_FLOAT, rank, 0, MPI_COMM_WORLD, &status);
    if (rank == 0) {
        cout << "B\n";
        for (int i = 0; i < m; i++) {

            cout << left
                << setw(5)
                << b[i]
                << " ";

        }
    }
    MPI_Type_free(&xpose);
    MPI_Finalize();
}