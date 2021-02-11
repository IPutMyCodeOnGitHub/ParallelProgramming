#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <pthread.h>
#include <math.h>
#include <iostream>
#include <fstream>

#define EPS 1e-15

using namespace std;

double fun(int k, int n, int i, int j);
int read_matrix_from_file(double* matrix, int n, string filename);
void generate_matrix_from_formula(double* matrix, int n, int k);
void print(double* matrix, int n, int m, ostream& of);
double mult_err(double* mat1, double* mat2, int n);
double norm(double* mat, int n);
int ind(int i, int j);
int get_inverse(double* mat, double* res, double* d, int n, int n_threads, int thread_i, int* res_code);

int read_matrix_from_file(double* matrix, int n, string filename) {
    ifstream in;
    in.open(filename);
    if (!in) {
        cout << "Error opening file" << endl;
        return -1;
    }
    double tmp;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (!(in >> tmp)) {
                cout << "Error reading file" << endl;
                return -1;
            }
        }
        for (int j = i; j < n; ++j) {
            if (!(in >> tmp)) {
                cout << "Error reading file" << endl;
                return -1;
            }
            matrix[ind(i, j)] = tmp;
        }
    }
    in.close();
    return 1;
}

void generate_matrix_from_formula(double* matrix, int n, int k) {
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            matrix[ind(i, j)] = fun(k, n, i, j);
        }
    }
}

void print(double* matrix, int n, int m, ostream& of) {
    for (int i = 0; i < min(n, m); ++i) {
        for (int j = 0; j < min(n, m); ++j) {
            of << matrix[ind(i, j)] << " ";
        }
        of << endl;
    }
}

double fun(int k, int n, int i, int j) {
    switch(k) {
        case 1: {
            return n - max(i, j) + 1;
        }
        case 2: {
            return max(i, j) + 1;
        }
        case 3: {
            return abs(i - j);
        }
        case 4: {
            return 1.0 / (i + j + 1);
        }
        default: {
            return 0;
        }
    }
}

long int get_full_time(void)
{
	struct timeval buf;
	gettimeofday(&buf, 0);
	return buf.tv_sec * 100 + buf.tv_usec/10000;
}
int ind(int i, int j) {
    if (i < j) return j * (j + 1) / 2 + i;
	else return i * (i + 1) / 2 + j;
}

double mult_err(double* mat1, double* mat2, int n) {
    double norm = 0;
    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double s = 0;
			for (int k = 0; k < n; ++k)
				s += mat1[ind(i, k)] * mat2[ind(k, j)];

            if (i == j)
                 s -= 1;

			norm += s * s;
		}
	}
	return sqrt(norm);
}

double norm(double* mat, int n) {
    double max = 0, sum = 0;
    for (int j = 0; j < n; ++j) {
        sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += abs(mat[ind(i, j)]);
        }
        if (sum > max) max = sum;
    }
    return max;
}

void synchronize(int total_threads) {
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}

/*
int ind(int i, int j) {
    if (i < j) return j * (j + 1) / 2 + i;
	else return i * (i + 1) / 2 + j;
}
*/

int get_inverse(double* mat, double* res, double* d, int n, int n_threads, int thread_i, int* res_code) {
    double s = 0;
	double no = 0;
	int row_from = 0, row_to = 0;

	if (thread_i == 0){
		s = 0;
		no = norm(mat, n);
		row_from = 0; row_to = 0;
	}

	row_from = n * thread_i / n_threads;
	row_to = n * (thread_i +  1) / n_threads;

	for (int i = row_from; i < row_to; ++i) {
		for (int j = i; j < n; ++j) {
			res[j * (j + 1) / 2 + i] = mat[j * (j + 1) / 2 + i];
		}
	}

	synchronize(n_threads);

	for (int k = 0; k < n; ++k) {
		if (thread_i == 0) {

			if (thread_i == 0 & fabs(res[k * (k + 1) / 2 + k] / no) < EPS) {
				*res_code = -1;
			}

			if (res[k * (k + 1) / 2 + k] < 0) {
				d[k] = -1;
				res[k * (k + 1) / 2 + k] = -res[k * (k + 1) / 2 + k];
			} else {
				d[k] = 1;
			}
			res[k * (k + 1) / 2 + k] = sqrt(res[k * (k + 1) / 2 + k]);

			for (int j = k + 1; j < n; ++j) {
				res[j * (j + 1) / 2 + k] /= res[k * (k + 1) / 2 + k] * d[k];
			}
		}


		synchronize(n_threads);
		if (*res_code == -1) return -1;

		row_from = (n - k - 1) * thread_i / n_threads + k + 1;
    	row_to = (n - k - 1) * (thread_i + 1) / n_threads + k + 1;

		for (int i = row_from; i < row_to; ++i) {
			for (int j = i; j < n; ++j) {
				res[j * (j + 1) / 2 + i] -= res[i * (i + 1) / 2 + k] * res[j * (j + 1) / 2 + k] * d[k];
			}
		}
		synchronize(n_threads);
	}

	//synchronize(n_threads);

	row_from = n * thread_i / n_threads;
    row_to = n * (thread_i + 1) / n_threads;

    for (int i = row_from; i < row_to; ++i) {

		for (int j = i; j >= 0; --j) {
			s = (double)(i == j);
			for (int k = j + 1; k <= i; ++k)
				s -= mat[i * (i + 1) / 2 + k] * res[k * (k + 1) / 2 + j];

			mat[i * (i + 1) / 2 + j] = s / res[j * (j + 1) / 2 + j];
		}
	}


	synchronize(n_threads);

	row_from = n * thread_i / n_threads;
    row_to = n * (thread_i + 1) / n_threads;
	for (int i = row_from; i < row_to; ++i) {
		for (int j = 0; j <= i; ++j)
		{
			s = 0.0;
			for (int k = i; k < n; ++k)
				s += d[k] * mat[k * (k + 1) / 2 + j] * mat[k * (k + 1) / 2 + i];
			res[i * (i + 1) / 2 + j] = s;
		}

		for (int j = i + 1; j < n; ++j)
		{
			s = 0.0;
			for (int k = j; k < n; ++k)
				s += d[k] * mat[k * (k + 1) / 2 + j] * mat[k * (k + 1) / 2 + i];
			res[j * (j + 1) / 2 + i] = s;
		}
	}
	synchronize(n_threads);
    return 1;
}

/*
int cholesky_decomp(double* mat, double* res, double* d, int n, int n_threads, int thread_i) {
    double s = 0;
    double no = norm(mat, n);
    for (int i = 0; i < n; ++i)
		for (int j = i; j < n; ++j)
		{
			s = mat[j * (j + 1) / 2 + i];
			for (int k = 0; k < i; ++k)
				s -= res[i * (i + 1) / 2 + k] * res[j * (j + 1) / 2 + k] * d[k];


			if (i == j)
			{
				if (s < 0)
				{
					s = -s;
					d[i] = -1.0;
				}
				else
					d[i] = 1.0;

				if (fabs(s / no) < EPS)
					return -1;

				s = sqrt(s);
				res[i * (i + 1) / 2 + i] = s;
			}
			else
				res[j * (j + 1) / 2 + i] = s / (res[i * (i + 1) / 2 + i] * d[i]);
		}
    return 1;
}*/

typedef struct
{
	int n;
	double* mat;
	double* d;
    double* res;
	int n_threads;
    int thread_i;
    int* code;
} ARGS;

void* solve(void* arg) {
    ARGS* args = (ARGS*) arg;
    get_inverse(args->mat, args->res, args->d, args->n, args->n_threads, args->thread_i, args->code);
	return NULL;
}

int main(int argc, char * argv[]) {
    string filename = "";
    int n, m, k, n_threads;
    if (argc == 6) {
        n = stoi(argv[1]);//size of matrix
        m = stoi(argv[2]);//size of part matrix (=n by default)
        k = stoi(argv[3]);//type of gen function
        n_threads = stoi(argv[4]);
        filename = argv[4];
    } else if (argc == 5) {
        n = stoi(argv[1]);
        m = stoi(argv[2]);
        k = stoi(argv[3]);
        n_threads = stoi(argv[4]);

    } else {
        cout << "Usage: ./main n m k filename or ./main n m k when k != 0";
        return 0;
    }

    double* matrix;
    double* d;
    double* res;
    pthread_t* threads;
    ARGS* args;
    try {
        res = new double[n * (n + 1) / 2];
        matrix = new double[n * (n + 1) / 2];
        d = new double[n];
        threads = (pthread_t*) malloc(n_threads * sizeof(pthread_t));
        args = (ARGS*) malloc(n_threads * sizeof(ARGS));
    } catch (const bad_alloc) {
        cout << "memory error" << endl;
        return 0;
    }

    if (k != 0) {
        generate_matrix_from_formula(matrix, n, k);
    } else {
        if (read_matrix_from_file(matrix, n, filename) != 1) {
            delete[] matrix; delete[] res; delete[] d;
            free(threads); free(args);
            return 0;
        }
    }

    cout << "MATRIX: " << endl;
    print(matrix, n, m, cout);
    cout << endl;

    int res_code = 0;
    for (int i = 0; i < n_threads; i++) {
		args[i].n = n;
		args[i].mat = matrix;
		args[i].d = d;
        args[i].res = res;
		args[i].thread_i = i;
		args[i].n_threads = n_threads;
        args[i].code = &res_code;
	}

    long int t_start = get_full_time();


    for (int i = 0; i < n_threads; i++)
		if (pthread_create(threads + i, 0, solve, args + i))
		{
			cout << "Error creating thread" << endl;
			delete[] matrix; delete[] res; delete[] d;
			if (threads) free(threads);
			if (args) free(args);
			return 0;
		}


	for (int i = 0; i < n_threads; i++)
		if (pthread_join(threads[i], 0))
		{
			cout << "Error joining thread" << endl;
			delete[] matrix; delete[] res; delete[] d;
            if (threads) free(threads);
			if (args) free(args);
			return 0;
		}

   //get_inverse(matrix, res, d, n, n_threads, 0);

    cout << "code: " << res_code << endl;
    if (res_code == -1) {
        cout << "det is zero" << endl;
         delete[] matrix; delete[] res;
         delete[] d;
         free(threads); free(args);
         return 0;
    }

    //get_inverse(matrix, res, d, n, n_threads, 0);
    cout << "THREADS: " << n_threads << endl;
    cout << "TIME: " << (double) (get_full_time() - t_start) << endl;

    cout << "A: " << endl;
    print(matrix, n, m, cout);
    cout << "RES: " << endl;
    print(res, n, m, cout);
    cout << endl;

    if (k != 0) {
        generate_matrix_from_formula(matrix, n, k);
    } else {
        if (read_matrix_from_file(matrix, n, filename) != 1) {
            delete[] matrix; delete[] res; delete[] d;
            free(threads); free(args);
            return 0;
        }
    }

    cout << "NORM: " << mult_err(matrix, res, n) << endl;
    delete[] matrix; delete[] res; delete[] d;
    free(threads); free(args);
}
