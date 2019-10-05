
#include "main_program.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>


void Ricker_amplitude_cpu(const float *F, float *force, int iter, float tao, float f, int n_size) {
	for (int i=0; i<n_size; ++i) {
		const float pi = 3.141592;
		float t = tao * iter - 1.0 / std::sqrt(2 * pi * pi * f * f); 
		float A = (1 - 2 * pi*pi * f*f * t*t) * exp(-pi*pi * f*f * t*t);

		force[i] =  F[i] * A;
	}
}

void make_A_x_y_cpu(const float *MasMatrix, float *StifMatrix,
					const float *ForceMatrix, float *U, float *V,
					int n_size, int step) {
	float tao = 0.001;
	clock_t start = clock();

	float *Alpha = new float[step*n_size];
	for (int i=0; i<n_size; ++i) {
		Alpha[i] = 0;
	}

	float *tmp_force = new float[n_size];
	int iter = 0;

	for (; iter<step-1; ++iter) {
		for (int i=0; i<n_size; ++i) {
			U[(iter+1)*n_size+i] = U[iter*n_size+i] + V[iter*n_size+i] * tao + Alpha[iter*n_size+i] * tao * tao / 2;

			float tmp = 0;
			for (int j=0; j<n_size; ++j) {
				if (StifMatrix[i*n_size+j] <= 0.000001 && StifMatrix[i*n_size+j] >= -0.000001) {
					StifMatrix[i*n_size+j] = 0;
				}
				tmp += StifMatrix[i*n_size+j] * U[(iter+1)*n_size+j]; // K[i*n_size+j]
			}
	
			Ricker_amplitude_cpu(ForceMatrix, tmp_force, iter, tao, 20, n_size);
			if (MasMatrix[i] == 0) {
				Alpha[(iter+1)*n_size+i] = 0;
			} else {
				Alpha[(iter+1)*n_size+i] = (tmp_force[i] - tmp) / MasMatrix[i];
			}

			V[(iter+1)*n_size+i] = V[iter*n_size+i] + tao * (Alpha[iter*n_size+i] + Alpha[(iter+1)*n_size+i]) / 2;
		}
	}

	clock_t end = clock();
	double seconds = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "CPU timing: (sec)" << seconds << std::endl;

	delete[] Alpha, tmp_force;
}


__device__ void Ricker_amplitude(float *F, float *force, int iter, float tao, float f, int offset) {
	const float pi = 3.141592;
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (offset != 0) {
		i += offset;
	}
	float t = tao * iter - 1.0 / std::sqrt(2 * pi * pi * f * f); 
	float A = (1 - 2 * pi*pi * f*f * t*t) * exp(-pi*pi * f*f * t*t);

	force[i] =  F[i] * A;
	//if (force[i] <= 0.000001 && force[i] >= -0.000001) {
	//	force[i] = 0;
	//}
}


__global__ void solving_system(float *M, float *K, float *F, float *U, float *V,
							   float *alpha, int iter, float tao, int n_size, float *tmp_force, int offset) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (offset != 0) {
		i += offset;
	}

	if (iter == 0) {
		alpha[i] = 0;
	}
	U[(iter+1)*n_size+i] = U[iter*n_size+i] + V[iter*n_size+i] * tao + alpha[iter*n_size+i] * tao * tao / 2;

	float tmp = 0;
	for (int j=0; j<n_size; ++j) {
		if (K[i*n_size+j] <= 0.000001 && K[i*n_size+j] >= -0.000001) {
			K[i*n_size+j] = 0;
		}
		tmp += K[i*n_size+j] * U[(iter+1)*n_size+j]; // K[i*n_size+j]
	}
	
	Ricker_amplitude(F, tmp_force, iter, tao, 20, offset);
	if (M[i] == 0) {
		alpha[(iter+1)*n_size+i] = 0;
	} else {
		alpha[(iter+1)*n_size+i] = (tmp_force[i] - tmp) / M[i];
	}
	//if (alpha[(iter+1)*n_size+i] <= 0.000001 && alpha[(iter+1)*n_size+i] >= -0.000001) {
	//	alpha[(iter+1)*n_size+i] = 0;
	//}

	V[(iter+1)*n_size+i] = V[iter*n_size+i] + tao * (alpha[iter*n_size+i] + alpha[(iter+1)*n_size+i]) / 2;
}


// also later need to recompute ForceMatrix  f = f(t)
// maybe make function which create all needed matrix except f
void make_A_x_y(const float *MasMatrix, const float *StifMatrix, const float *ForceMatrix,
				float *U, float *V, int n_size, int step) {
	float tao = 0.001;
	float theta1 = 1, theta2 = 0;
	int size = n_size * n_size * sizeof(float);

	float *d_M, *d_K;//, *d_F;//, *h_A;

	cudaEvent_t start;
    cudaEvent_t stop;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start, 0);
	// cudaMalloc((void**)&d_M, size);
	cudaMalloc((void**)&d_M, size / n_size);
	cudaMalloc((void**)&d_K, size);
	//cudaMalloc((void**)&d_A, size);

	// want to get A matrix in A*x = y:
	// cudaMemcpy(d_M, MasMatrix, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_M, MasMatrix, size / n_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_K, StifMatrix, size, cudaMemcpyHostToDevice);

	int block_size = 637;
	int grid_size = 26; // (64 * 64 + block_size - 1) / block_size;
	
	float *d_F, *d_u, *d_v, *d_alpha;
	float *d_tmp_force;
	cudaMalloc((void**)&d_tmp_force, size / n_size);
	cudaMalloc((void**)&d_F, size / n_size);
	cudaMalloc((void**)&d_u, size / n_size * step);
	cudaMalloc((void**)&d_v, size / n_size * step);
	cudaMalloc((void**)&d_alpha, size / n_size * step);

	cudaMemcpy(d_u, U, size / n_size * step, cudaMemcpyHostToDevice);
	cudaMemcpy(d_v, V, size / n_size * step, cudaMemcpyHostToDevice);
	cudaMemcpy(d_F, ForceMatrix, size / n_size, cudaMemcpyHostToDevice);
	
	// actually compute Y here:
	// this part need to be iterated, iter is num of iteration for
	// choosing right column in U and V
	int iter = 0;

	// CHANGE NUM OF THREADS ITS IPMORTANT
	block_size = 289;
	grid_size = 8; //(n_size / 256) + 1;
	// cudaEventRecord(start, 0);
	// also need to evaluate velocity
	float *Alpha = new float[step*n_size];
	for (; iter<step-1; ++iter) {
		solving_system<<<grid_size, block_size>>>(d_M, d_K, d_F, d_u, d_v, d_alpha, iter,
												  tao, n_size, d_tmp_force, 0);
		//solving_system<<<grid_size, block_size>>>(d_M, d_K, d_F, d_u, d_v, d_alpha, iter,
		//										  tao, n_size, d_tmp_force, 9216);
	}
	// cudaEventRecord(stop, 0);

	cudaMemcpy(U, d_u, size / n_size * step, cudaMemcpyDeviceToHost);
	cudaMemcpy(V, d_v, size / n_size * step, cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);
	float time = 0;
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    std::cout << "GPU compute time: " << time << std::endl;

	// dont forget to free all memory !!!!
	cudaFree(d_u);
	cudaFree(d_v);
	cudaFree(d_M);
	cudaFree(d_K);
	cudaFree(d_F);
	cudaFree(d_alpha);
	cudaFree(d_tmp_force);
}

int main()
{
	// change force in main_program: 50 line
	// change mesh in load_mesh: 54 line
	// change num of threads: 241 line
	// change filename to result
	// dont forget to delete first line in elem.txt and node.txt
	int n_size = 1156 * 2;
	int step = 800;

	float *MasMatrix = new float[n_size];
	float *StifMatrix = new float[n_size*n_size];
	float *ForceMatrix = new float[n_size];
	main_program(MasMatrix, StifMatrix, ForceMatrix, n_size);

	// boarder condition that on y == -2 and y == 2: displacment == 0
	float y_bord = 2.5, x_bord = 5;
	for (int i=0; i<Nodes.size(); ++i) {
		if (/*Nodes[i].y == y_bord || */Nodes[i].y == -y_bord
			|| Nodes[i].x == x_bord || Nodes[i].x == -x_bord) {
		// if (Nodes[i].x == 2 || Nodes[i].x == -2) {
			int k = 2 * (Nodes[i].nid - 1);
			ForceMatrix[k] = 0;
			MasMatrix[k] = 0;
			for (int j=0; j<n_size; ++j) {
				StifMatrix[k*n_size+j] = 0;
				StifMatrix[j*n_size+k] = 0;
				if (k == j) {
					StifMatrix[k*n_size+j] = 1;
				}
			}
		}

		if (/*Nodes[i].y == y_bord || */Nodes[i].y == -y_bord
			|| Nodes[i].x == x_bord || Nodes[i].x == -x_bord) {
			int k = 2 * (Nodes[i].nid);
			ForceMatrix[k] = 0;
			MasMatrix[k] = 0;
			for (int j=0; j<n_size; ++j) {
				StifMatrix[k*n_size+j] = 0;
				StifMatrix[j*n_size+k] = 0;
				if (k == j) {
					StifMatrix[k*n_size+j] = 1;
				}
			}
		}
	}

	//Elements.clear();
	Nodes.clear();
	std::cout << std::endl;

	float *U = new float[step*n_size];
	float *V = new float[step*n_size];

	for (int i=0; i<step; ++i) {
		for (int j=0; j<n_size; ++j) {
			U[i*n_size+j] = 0;
			V[i*n_size+j] = 0;
			if (i == 0) {
				// need to write U0 and V0
				U[j] = 0;
				V[j] = 0;
			}
		}
	}

	make_A_x_y(MasMatrix, StifMatrix, ForceMatrix, U, V, n_size, step);

	make_A_x_y_cpu(MasMatrix,StifMatrix, ForceMatrix, U, V, n_size, step);

	/*std::ofstream out;
	out.open("result_velocities_XxX_xxx_nu20_f4000_newriecker_mu0_step800_50try_bordercond.txt");
	for (int i=0; i<n_size; ++i) {
		out << (int)(i / 2) << " ";
		for (int j=0; j<step; ++j) {
			out << V[j*n_size+i] << " ";
		}
		out << '\n';
	}*/

	delete[] MasMatrix, StifMatrix, ForceMatrix;
	delete[] U, V;

	return 0;
}