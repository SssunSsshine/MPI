#include <iostream>
#include <mpi.h>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <functional>
/*20.	 аждый процесс обмениваетс€ сообщени€ми со всеми остальными и выводит номера тех процессов, 
у которых ни одно значение не превышает его минимальное значение.*/

using namespace std;


const int m = 5;
int TAG = 0;

int get_min(int  arr[m]) {
	int minimum = arr[0];
	for (int i = 1; i < m; i++)
		minimum = min(minimum, arr[i]);
	return minimum;
}

bool is_less(int min, int  buffer[m])
{
	bool is_more = true;	
	for (int i = 0; i < m && is_more; i++) {
		is_more = min >= buffer[i];
	}
	return is_more;
}

void print_array(int  message[m])
{
	for (int i = 0;i < m;i++)
		cout << left
		<< setw(10)
		<< message[i] << ' ';
	cout << "\n";
}

int main(int argc, char* argv[])
{
	std::vector<int> vmain;
	int procnum;
	int message[m];
	int myrank;
	MPI_Status status;

	bool isfirst = true;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);

	if (myrank % 2 == 0) {
		for (int i = 0; i < m; i++) {
			message[i] = myrank + (i + 1) * 10;
		}
	}
	else {
		for (int i = 0; i < m; i++) {
			message[i] = myrank + i;
		}
	}
	
	cout << myrank << " process: ";
	print_array(message);

	for (int i = 0; i < procnum; i++) {
		if (i != myrank) {
			MPI_Send(&message, m, MPI_INT, i, TAG, MPI_COMM_WORLD);
		}
	}

	cout << myrank << " process: ";

	int min = get_min(message);
	cout << "min " << min << " processes: ";
	for (int i = 0; i < procnum; i++) {
		if (i != myrank) {
			int buffer[m];
			MPI_Recv(&buffer, m, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);

			if (is_less(min, buffer)) {
				cout << i << " ";
			}
		}
	}
	cout << "\n\n";

	MPI_Finalize();
}


