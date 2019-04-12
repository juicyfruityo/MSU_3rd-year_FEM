
#include "main_program.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

#include <iostream>


__global__ void make_A_matrix(float *M, float *K, float *A, float theta2, float tao) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	A[i] = M[i] * 2.0/theta2 + tao * tao * K[i];
}

__global__ void make_Y_matrix(float *M, float *U, float *V, float *F, float *Y, float theta1,
							   float theta2, float tao, int iter) {
	int tx = threadIdx.x + blockIdx.x * blockDim.x;
	Y[tx] = 2/theta2 * (U[iter*64+tx] + theta1 * tao * V[iter*64+tx]);
	// need to use shared memory
	float res = 0;
	for (int i=0; i<64; ++i) {
		res += Y[i] * M[tx*64+i];
	}
	Y[tx] = res + F[tx];
}

// create permutation matrix for solving system
void createPermutationMatrix(int *h_PivotArray, float *P, int N) {
    int temp;
    // --- Create permutation matrix
	// float *P = new float[N * N * sizeof(float)];
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++)
			if (i == j) P[i * N + j] = 1;
			else P[i*N+j] = 0;
	}
	for (int j=0; j<N; j++) 
		for (int i=0; i<N; i++) {
			temp = P[i + j * N];
			P[i + j * N] = P[(h_PivotArray[i] - 1) + j * N];
			P[(h_PivotArray[i] - 1) + j * N] = temp;
		}
}

void solve_system(float *h_A, float *h_B, float *U, int iter) {
	const unsigned int N = 64; 
	cublasHandle_t handle;
	cublasCreate(&handle);
	/***********************/
	/* SETTING THE PROBLEM */
	/***********************/
	// --- Matrices to be inverted (only one in this example)
	// need to transpose this matrix
	for (int i=0; i<N; ++i) {
		for (int j=0; j<N; ++j) {
			float tmp = h_A[i*N+j];
			h_A[i*N+j] = h_A[j*N+i];
			h_A[j*N+i] = tmp;
		}
	}
	// h_X will be Un+1 in big matrix U
	float *h_X = new float[N];

	float *d_A;	cudaMalloc((void**)&d_A, N*N*sizeof(float));
	// d_B == h_Y in my code, just because of legacy (i stole this code)
	float *d_B;	cudaMalloc((void**)&d_B, N*sizeof(float));
	// d_X is Un+1
	// float *d_X;	cudaMalloc((void**)&d_X, 10*N*sizeof(float)); // 10 because of i got 10 iteration

	// --- Move the relevant matrices from host to device
	cudaMemcpy(d_A, h_A, N*N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_B, h_B, N*sizeof(float), cudaMemcpyHostToDevice);
	/********************/
	/* LU DECOMPOSITION */
	/********************/
	// --- Creating the array of pointers needed as input/output to the batched getrf
	float **h_inout_pointers = (float **)malloc(sizeof(float *));
	h_inout_pointers[0] = (float *)((char*)d_A);
 
	float **d_inout_pointers;
	cudaMalloc((void**)&d_inout_pointers, sizeof(float *));
	cudaMemcpy(d_inout_pointers, h_inout_pointers, sizeof(float *), cudaMemcpyHostToDevice);
	free(h_inout_pointers);

	int *d_PivotArray; cudaMalloc((void**)&d_PivotArray, N*sizeof(int));
	int *d_InfoArray;  cudaMalloc((void**)&d_InfoArray,  sizeof(int));
	
	int *h_PivotArray = (int *)malloc(N*sizeof(int));
	int *h_InfoArray  = (int *)malloc(sizeof(int));
	// actually making LU factorization
	cublasSgetrfBatched(handle, N, d_inout_pointers, N, d_PivotArray, d_InfoArray, 1);
	
    cudaMemcpy(h_InfoArray, d_InfoArray, sizeof(int), cudaMemcpyDeviceToHost);

	if (h_InfoArray[0]  != 0) {
		fprintf(stderr, "Factorization of matrix %d Failed: Matrix may be singular\n", h_InfoArray[0]);
		cudaDeviceReset();
		// here i need throw exception or something 
	}

	cudaMemcpy(h_A, d_A, N*N*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_PivotArray, d_PivotArray, N*sizeof(int),cudaMemcpyDeviceToHost);

	// --- The output factored matrix in column-major format
	// for (int i=0; i<N*N; i++) printf("A[%i]=%f\n", i, h_A[i]);

	// printf("\n\n");
	// --- The pivot array
	// for (int i=0; i<N; i++) printf("IPIV[%i]=%i\n", i, h_PivotArray[i]);

	float *P = new float[N * N * sizeof(float)];
	createPermutationMatrix(h_PivotArray, P, N);

	float *d_P; cudaMalloc((void**)&d_P, N * N * sizeof(float));

	// printf("\n\n");
	// --- The permutation matrix
	// for (int i=0; i<N; i++) 
		// for (int j=0; j<N; j++)
			// printf("P[%i, %i]=%f\n", i, j, P[j * N + i]);

	cudaMemcpy(d_P, P, N * N * sizeof(float), cudaMemcpyHostToDevice);
	
	// --- Now P*A=L*U
	//     Linear system A*x=y => P.'*L*U*x=y => L*U*x=P*y
	float alpha1  = 1.f;
	float beta1   = 0.f;
	// actually solving problem from here
	cublasSgemv(handle, CUBLAS_OP_N, N, N, &alpha1, d_P, N, d_B, 1, &beta1, d_B, 1);

	// cudaMemcpy(h_B,d_B,N*sizeof(float),cudaMemcpyDeviceToHost);
	
	// --- The result of P*y
	// printf("\n\n");
	// for (int i=0; i<N; i++) printf("(P*y)[%i]=%f\n", i, h_B[i]);

	const float alpha  = 1.f;

	// --- Function solves the triangulatr linear system with multiple right hand sides, function overrides b as a result 

	// --- Lower triangular part
	cublasStrsm(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,
				 CUBLAS_DIAG_UNIT, N, 1, &alpha, d_A, N, d_B, N);

	// --- Upper triangular part
	cublasStrsm(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N,
				CUBLAS_DIAG_NON_UNIT, N, 1, &alpha, d_A, N, d_B, N);

	cudaMemcpy(h_X,d_B,N*sizeof(float),cudaMemcpyDeviceToHost);
	
	for (int i=0; i<N; ++i) {
		U[(iter+1)*N+i] = h_X[i];
	}

	delete[] h_X, P;
	free(d_PivotArray);
	free(d_InfoArray);
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_inout_pointers);
	cudaFree(d_PivotArray);
	cudaFree(d_InfoArray);
	cudaFree(d_P);
	// --- The output inverted matrix in column-major format
	// printf("\n\n");
	// for (int i=0; i<N; i++) printf("B[%i]=%f\n", i, h_B[i]);

}


// also later need to recompute ForceMatrix  f = f(t)
// maybe make function which create all needed matrix except f
void make_A_x_y(const float *MasMatrix, const float *StifMatrix, const float *ForceMatrix,
				float *U, float *V) {
	float tao = 0.001;
	float theta1 = 1, theta2 = 1;
	int size = 64 * 64 * sizeof(float);

	float *d_m, *d_k, *d_A, *h_A;
	h_A = new float[64*64];

	cudaMalloc((void**)&d_m, size);
	cudaMalloc((void**)&d_k, size);
	cudaMalloc((void**)&d_A, size);

	// want to get A matrix in A*x = y:
	cudaMemcpy(d_m, MasMatrix, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_k, StifMatrix, size, cudaMemcpyHostToDevice);

	int block_size = 512;
	int grid_size = 8; // (64 * 64 + block_size - 1) / block_size;
	make_A_matrix<<<grid_size, block_size>>>(d_m, d_k, d_A, theta2, tao);
	
	cudaFree(d_k);

	// where Y: A * x = Y
	float *d_f, *d_u, *d_v, *d_Y, *h_Y;
	h_Y = new float[64];
	cudaMalloc((void**)&d_f, size / 64);
	cudaMalloc((void**)&d_u, size / 64 * 10);
	cudaMalloc((void**)&d_v, size / 64 * 10);
	cudaMalloc((void**)&d_Y, size / 64);

	cudaMemcpy(d_u, U, size / 64, cudaMemcpyHostToDevice);
	cudaMemcpy(d_v, V, size / 64, cudaMemcpyHostToDevice);
	cudaMemcpy(d_f,ForceMatrix, size / 64, cudaMemcpyHostToDevice);
	
	// actually compute Y here:
	// this part need to be iterated, iter is num of iteration for
	// choosing right column in U and V
	int iter = 0;
	block_size = 64;
	grid_size = 1;
	make_Y_matrix<<<grid_size, block_size>>>(d_m, d_u, d_v, d_f, d_Y, theta1, theta2, tao, iter);

	// next need to solve A * x = Y (x == Un, A == d_A, Y == d_Y)
	cudaMemcpy(h_A, d_A, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_Y, d_Y, size/64, cudaMemcpyDeviceToHost);

	solve_system(h_A, h_Y, U, iter);
	// also need to evaluate velocity

	// dont forget to free all memory !!!!
}

int main()
{
	float *MasMatrix = new float[64*64];
	float *StifMatrix = new float[64*64];
	float *ForceMatrix = new float[64];
	main_program(MasMatrix, StifMatrix, ForceMatrix);

	float *U = new float[10*64];
	float *V = new float[10*64];

	for (int i=0; i<10; ++i) {
		// U[i] = new float[64];
		// V[i] = new float[64];
		for (int j=0; j<64; ++j) {
			U[i*64+j] = 0;
			V[i*64+j] = 0;
			if (i == 0) {
				// need to write U0 and V0
				U[j] = j;
				V[j] = j;
			}
		}
	}

	make_A_x_y(MasMatrix, StifMatrix, ForceMatrix, U, V);

	return 0;
}