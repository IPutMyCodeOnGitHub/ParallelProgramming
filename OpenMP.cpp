#include <iostream>
#include <bits/stdc++.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <omp.h>

#define N 10//N of matrix
//int N=0;

bool DEBUG=false;

const int Tag = 0;
const int root = 0;

using namespace std;



// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
void getCofactor(int A[N][N], int temp[N][N], int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			// Copying into temporary matrix only those element
			// which are not in given row and column
			if (row != p && col != q)
			{
				temp[i][j++] = A[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
n is current dimension of A[][]. */
int determinant(int A[N][N], int n)
{
	int D = 0; // Initialize result

	// Base case : if matrix contains single element
	if (n == 1)
		return A[0][0];

	int temp[N][N]; // To store cofactors

	int sign = 1; // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of A[0][f]
		getCofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(int A[N][N],int adj[N][N])
{
	if (N == 1)
	{
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1, temp[N][N];

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			// Get cofactor of A[i][j]
			getCofactor(A, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i+j)%2==0)? 1: -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign)*(determinant(temp, N-1));
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(int A[N][N], float inverse[N][N])
{
	// Find determinant of A[][]
	int det = determinant(A, N);
	if (det == 0)
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	int adj[N][N];
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			inverse[i][j] = adj[i][j]/float(det);

	return true;
}

// Generic function to display the matrix. We use it to display
// both adjoin and inverse. adjoin is integer matrix and inverse
// is a float.
template<class T>
void display(T A[N][N])
{
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}
}

// Driver program
int main(int argc,char *argv[])
{



  int A[N][N] = { {5, -2, 2, 7},
                {1, 0, 0, 3},
                {-3, 1, 5, 0},
                {3, -1, -9, 4}};


  if(argc!=2){
    cout << "Not enough arguments. Must be 2. Restart program" << endl;
      exit(1);
    }
  int n=stoi(argv[1]);
  if (N%n!=0 || n>N){
    cout << "Quantity of threads must be one of dividers of "<< N <<". Restart program with one of these arguments:" << endl;
    for (int i = 1; i <= N; i++){
       if (N % i == 0)
           cout << i << endl;
         }
    exit(0);
  }
  omp_set_num_threads(n);
  if(!DEBUG){
    for (int i = 0; i<N; i++){
      for (int j = 0; j<N; j++){
        A[i][j] = rand();
      }
    }
  }

  int adj[N][N]; // To store adjoint of A[][]
  float inv[N][N]; // To store inverse of A[][]

  auto start=chrono::system_clock::now();
  
#pragma omp parallel shared(adj)
{
  adjoint(A, adj);
}

#pragma omp parallel shared(inv)
{
  inverse(A, inv);
}
  auto end = chrono::system_clock::now();
  chrono::duration<double> elapsed = end - start;

	cout << "Input matrix is :" << endl;
	display(A);
	cout << "\nThe Adjoint is :" << endl;
	display(adj);
	cout << "\nThe Inverse is :" << endl;
	if (inv)
		display(inv);
  else
    cout << "Error has been occured" << endl;
  cout << "Elapsed time: " << elapsed.count() << "s" << endl;



	return 0;
}
