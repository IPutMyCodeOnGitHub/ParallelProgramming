import numpy as np
from numpy.linalg import inv
from scipy.linalg import det
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

comm.Barrier()
t_start = MPI.Wtime()

n=None

if rank==0:
    n = int(input("Enter size of matrix (n*n)\n"))

n = comm.bcast(n, root=0)

if rank==0:
    matrix = np.zeros((n,n))
    #print(matrix)
    determinant = np.zeros((2))
    for i in range (0,n):
        for j in range(0,n):
            #print("Enter the value for the field: n = ", i+1, ", n = ", j+1,"\n")
            #matrix[i][j] = float(input())
            matrix[i][j]=np.random.sample()

    comm.Send(matrix, dest=1, tag=0)
    comm.Send(matrix, dest=2, tag=0)
    comm.Recv(determinant, source=1, tag=0)
    comm.Recv(matrix, source=2, tag=0)
    if(determinant[0]==0):
        print("There is no inverse matrix since the determinant is 0!")
    else:
        print("Inverse matrix:\n",matrix)

if rank==1:
    matrix = np.zeros((n,n))
    determinant = np.zeros((2))
    comm.Recv(matrix, source=0, tag=0)
    determinant = np.array([[det(matrix)],[0]])
    comm.Send(determinant, dest=0)
if rank==2:
    matrix = np.zeros((n,n))
    comm.Recv(matrix, source=0, tag=0)
    matrix = inv(matrix)
    comm.Send(matrix, dest=0)
#else:
#    exit()
comm.Barrier()
t_diff = MPI.Wtime() - t_start

if comm.rank==0:
    print ("Elapsed time "+str(t_diff)+"s");
