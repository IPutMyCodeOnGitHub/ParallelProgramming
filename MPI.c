#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define MAX_ARRAY_SIZE 10
#define TRUE 1
#define FALSE 0

int FindDetPara(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int my_rank, int p);
int FindDet(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE] , int SizeOfMatrix);
int Find2x2(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE] , int SizeOfMatrix);
int FindMinorMatrix(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int xpos, int ypos, int MinorMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE]);
void FindMatrixOfMinors(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int MatrixOfMinors[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int my_rank, int p);
void FindTransposedMatrix(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int TransposedMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE]);
void ChangeSigns(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix);
void PrintDecimalMatrix(float Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix);
void PrintMatrix(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix);

void FindInverseMatrix(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, float InverseMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int my_rank, int p){
/*	int MatrixOfMinors[SizeOfMatrix - 1][SizeOfMatrix - 1];
	int TransposedMatrix[SizeOfMatrix - 1][SizeOfMatrix - 1]; */
	int MatrixOfMinors[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
	int TransposedMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
	int DetMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
	int Det, y , x, pos, i, power, source, tag, tmpDet, dest;
	MPI_Status status;

	/*each process calculates their part of det*/
	Det = FindDetPara(Matrix, SizeOfMatrix, my_rank, p);
	FindMatrixOfMinors(Matrix, SizeOfMatrix, MatrixOfMinors, my_rank, p);

	/*Det = Find2x2(Matrix, SizeOfMatrix);*/
	if(my_rank == 0){

	if (SizeOfMatrix>2) {
		/*FindMatrixOfMinors(Matrix, SizeOfMatrix, MatrixOfMinors);
		printf("Matrix of minors:\n");
		PrintMatrix(MatrixOfMinors, SizeOfMatrix);*/
		ChangeSigns(MatrixOfMinors, SizeOfMatrix);
		//printf("Matrix of minors sign changed:\n");
		//PrintMatrix(MatrixOfMinors, SizeOfMatrix);
		FindTransposedMatrix(MatrixOfMinors, SizeOfMatrix, TransposedMatrix);
		//printf("Matrix transposed:\n");
		//PrintMatrix(TransposedMatrix, SizeOfMatrix);


		for(y = 0; y < SizeOfMatrix; y++){
			for(x = 0; x < SizeOfMatrix; x++){
				/*printf("Transpose=%d, Det=%d\n", TransposedMatrix[x][y], Det);*/
				InverseMatrix[x][y] = (float)TransposedMatrix[x][y] * (1 / (float)Det);
			}
		}
	} else {
		/* special case code for 2x2 matrix should go here */
	}

	}
}

int FindDetPara(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int my_rank, int p){
	int DetMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
	int Det, pos, i, power, source, tag, tmpDet, dest;
	MPI_Status status;

	Det = 0;
	tmpDet = 0;
	for(i = 0; ((my_rank + 1) + (p*i)) <= SizeOfMatrix; i++){
		pos = ((my_rank + 1) + (p*i)) - 1;
		if(pos % 2 == 0 | pos == 0){
			power = 1;
		}else{
			power = -1;
		}
		FindMinorMatrix(Matrix, SizeOfMatrix, pos, 0, DetMatrix);
		Det += power * Matrix[pos][0] * Find2x2(DetMatrix, SizeOfMatrix - 1);
	}
	/*Collect the Dets */
	if(my_rank == 0){
		for(source = 1; source < p; source++){
			tag = 2;
			MPI_Recv(&tmpDet, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			//printf("Recieve Det %d from process %d\n", tmpDet, source);
			Det += tmpDet;
		}
	//printf("Det = %d\n" , Det);
	}else{
		dest = 0;
		tag = 2;
		MPI_Send(&Det, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
	}

	return Det;
}

int FindDet(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix){
/*	int DetMatrix[SizeOfMatrix - 2][SizeOfMatrix - 2]; */
	int DetMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
	int i, Det=0;
	int power = -1;
	for(i = 0; i < SizeOfMatrix; i++){
                power *= -1;
		FindMinorMatrix(Matrix, SizeOfMatrix, i, 0, DetMatrix);
		/*PrintMatrix(DetMatrix, SizeOfMatrix - 1);*/
		Det += power * Matrix[i][0] * Find2x2(DetMatrix, SizeOfMatrix - 1);
	}
	return Det;

}

int Find2x2(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE] , int SizeOfMatrix){
	int Det=0;
	if(SizeOfMatrix == 2){
		Det = (Matrix[0][0] * Matrix[1][1]) - (Matrix[0][1] * Matrix[1][0]);
	}else{
		Det = FindDet(Matrix, SizeOfMatrix);
	}
	return Det;
}

int FindMinorMatrix(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int xpos, int ypos, int MinorMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE]){
	int relx=0, rely=0, y, x;
/*	int MinorMatrix[SizeOfMatrix - 2][SizeOfMatrix - 2]; */
/*	int MinorMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE]; */
	for(y = 0; y < SizeOfMatrix; y++){
		if(y != ypos){
			for(x = 0; x < SizeOfMatrix; x++){
				if(x != xpos){
					MinorMatrix[relx][rely] = Matrix[x][y];
					relx += 1;
				}
			}
			relx = 0;
			rely += 1;
		}
	}
}
void FindMatrixOfMinors(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int MatrixOfMinors[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int my_rank, int p){
	int x , y;
	int MinorMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];

	for(y = 0; y < SizeOfMatrix; y++){
		for(x = 0; x < SizeOfMatrix; x++){
			FindMinorMatrix(Matrix, SizeOfMatrix, x, y, MinorMatrix);
			MatrixOfMinors[x][y] = FindDetPara(MinorMatrix, SizeOfMatrix - 1, my_rank, p);
		}
	}
}

void FindTransposedMatrix(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix, int TransposedMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE]){
	int x , y;
	for(y = 0; y < SizeOfMatrix; y++){
		for(x = 0; x < SizeOfMatrix; x++){
			TransposedMatrix[x][y] = Matrix[y][x];
		}
	}
}

void ChangeSigns(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix){
	int i = 1, row = 1;
	int x , y;
	for(y = 0; y < SizeOfMatrix; y++){
		for(x = 0; x < SizeOfMatrix; x++){
			if(i == TRUE){
				i = FALSE;
			}else{
				Matrix[x][y] = Matrix[x][y] * -1;
				i = TRUE;
			}
		}
		if(row == TRUE){
			row = FALSE;
			i = FALSE;
		}else{
			row = TRUE;
			i = TRUE;
		}
	}
}

void PrintDecimalMatrix(float Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix) {
	int x, y;
	for(y = 0; y < SizeOfMatrix; y++){
		for(x = 0; x < SizeOfMatrix; x++){
			printf("%f ", Matrix[x][y]);
		}
		printf("\n");
	}

}

void PrintMatrix(int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE], int SizeOfMatrix) {
	int x, y;
	printf("SizeOfMatrix =%d\n", SizeOfMatrix);
	for(y = 0; y < SizeOfMatrix; y++){
		for(x = 0; x < SizeOfMatrix; x++){
			printf("%d ", Matrix[x][y]);
		}
		printf("\n");
	}

}

int main(int argc, char** argv){
	int SizeOfMatrix, y, x;
	int Matrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
	float InverseMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];

	/*MPI setup*/
	int my_rank, source, dest, tag, p;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);


	clock_t begin = clock();
	if(my_rank == 0){
		srand(time(NULL));
		int SizeOfMatrix=20;
		//printf("Enter size of matrix:\n");
		//scanf("%d" , &SizeOfMatrix);
		//printf("SizeOfMatrix = %d\n", SizeOfMatrix);
	/*	int Matrix[SizeOfMatrix - 1][SizeOfMatrix - 1]; */
	/*	float InverseMatrix[SizeOfMatrix - 1][SizeOfMatrix - 1]; */
		for(y = 0; y < SizeOfMatrix; y++){
			for(x = 0; x < SizeOfMatrix; x++){
				/*scanf("%d" , Matrix[x][y]);
				printf(",");*/
				Matrix[x][y] = rand() % 10;
			}
		}
	/*Send out matrix to processes */
	for(dest = 1; dest < p; dest++){
		tag = 0;
		MPI_Send(&Matrix,MAX_ARRAY_SIZE*MAX_ARRAY_SIZE,MPI_INT,dest,tag,MPI_COMM_WORLD);
		tag = 1;
		MPI_Send(&SizeOfMatrix,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
	}

	//printf("Calling print matrix\n");
  //PrintMatrix(Matrix, SizeOfMatrix);

	printf("Calling inverse matrix\n");
	FindInverseMatrix(Matrix, SizeOfMatrix, InverseMatrix, my_rank, p);

	printf("Calling print dec matrix\n");
  PrintDecimalMatrix(InverseMatrix, SizeOfMatrix);

	}
	else{
	source = 0;
	tag = 0;
	MPI_Recv(&Matrix,MAX_ARRAY_SIZE*MAX_ARRAY_SIZE,MPI_INT,source,tag,MPI_COMM_WORLD,&status);
	tag = 1;
	MPI_Recv(&SizeOfMatrix,1,MPI_INT,source,tag,MPI_COMM_WORLD,&status);
	FindInverseMatrix(Matrix, SizeOfMatrix, InverseMatrix, my_rank, p);
	}
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	if(my_rank == 0){
		printf("Elapsed time: %f\n",time_spent);
	}

	MPI_Finalize();
}
