#ifndef MAKE_LOCAL_MATRIX
#define MAKE_LOCAL_MATRIX

#include "load_mesh.h"
#include <cmath>


float basis_function(int num, float xi, float eta) {
  float result;

  switch (num) {
      case 0: result = (1-xi) * (1-eta) / 4;
              break;
      case 1: result = (1+xi) * (1-eta) / 4;
              break;
      case 2: result = (1+xi) * (1+eta) / 4;
              break;
      case 3: result = (1-xi) * (1+eta) / 4;
              break;
  }

  return result;
}

float Jacobian(element &elem, float xi, float eta) {
  float result=0;

  float dxdxi = (-1.0/2) * (1-eta)/2 * elem._node[0].x
              + 1.0/2 * (1-eta)/2 * elem._node[1].x
              + 1.0/2 * (1+eta)/2 * elem._node[2].x
              + (-1.0/2) * (1+eta)/2 * elem._node[3].x;

  float dxdeta = (-1.0/2) * (1-xi)/2 * elem._node[0].x
                + (-1.0/2) * (1+xi)/2 * elem._node[1].x
                + 1.0/2 * (1+xi)/2 * elem._node[2].x
                + 1.0/2 * (1-xi)/2 * elem._node[3].x;

  float dydxi = (-1.0/2) * (1-eta)/2 * elem._node[0].y
              + 1.0/2 * (1-eta)/2 * elem._node[1].y
              + 1.0/2 * (1+eta)/2 * elem._node[2].y
              + (-1.0/2) * (1+eta)/2 * elem._node[3].y;

  float dydeta = (-1.0/2) * (1-xi)/2 * elem._node[0].y
                + (-1.0/2) * (1+xi)/2 * elem._node[1].y
                + 1.0/2 * (1+xi)/2 * elem._node[2].y
                + 1.0/2 * (1-xi)/2 * elem._node[3].y;

  result = dxdxi * dydeta - dxdeta * dydxi;

  return result;
}

void mass_matrix_local(element &elem, float *Mloc) {
  float ro = 10000;
  std::vector<float> quad;
  quad.push_back(-1); quad.push_back(1);
  std::vector<float> quad_w;
  quad_w.push_back(1); quad_w.push_back(1);

  for (int i=0; i<8; ++i) {
      for (int j=0; j<8; ++j) {
          float tmp=0;
          // квадратуры гауса
          for (int k=0; k<quad.size(); ++k) {
              for (int l=0; l<quad.size(); ++l) {
                  tmp += basis_function(i/2, quad[k], quad[l])
                      * basis_function(j/2, quad[k], quad[l])
                      * Jacobian(elem, quad[k], quad[l])
                      * quad_w[k] * quad_w[l];
              }
          }
          Mloc[i*8+j] = tmp * ro;
      }
  }
}

void assembly_one_matrix(element &elem, float *Matrix, float *Mlocal, int n_size) {
  for (int i=1; i<=8; ++i) {
      int I = elem.num[(int)((i-1)/2) + 0]*2 - i%2;
      for (int j=1; j<=8; ++j) {
          int J = elem.num[(int)((j-1)/2) + 0]*2 - j%2;
          Matrix[(I-1)*n_size+(J-1)] += Mlocal[(i-1)*8+(j-1)];
      }
  }
}

void assembly_mass_matrix(element &elem, float *Matrix, float *Mlocal, int n_size) {
  for (int i=1; i<=8; ++i) {
      int I = elem.num[(int)((i-1)/2) + 0]*2 - i%2;
      for (int j=1; j<=8; ++j) {
          int J = elem.num[(int)((j-1)/2) + 0]*2 - j%2;
		  if (I == J) {
				Matrix[(I-1)] += Mlocal[(i-1)*8+(j-1)];
		  }
      }
  }
}

void assembly_force_matrix(element &elem, float *Matrix, float *Mlocal) {
  for (int i=1; i<=8; ++i) {
      int I = elem.num[(int)((i-1)/2) + 0]*2 - i%2;
      Matrix[(I-1)] += Mlocal[(i-1)];
  }
}

// Делаю B и  B_t
void make_grad_matrix(element &elem, float *B, float *B_t, float xi, float eta) {
  float dxdxi = (-1.0/2) * (1-eta)/2 * elem._node[0].x
              + 1.0/2 * (1-eta)/2 * elem._node[1].x
              + 1.0/2 * (1+eta)/2 * elem._node[2].x
              + (-1.0/2) * (1+eta)/2 * elem._node[3].x;

  float dxdeta = (-1.0/2) * (1-xi)/2 * elem._node[0].x
                + (-1.0/2) * (1+xi)/2 * elem._node[1].x
                + 1.0/2 * (1+xi)/2 * elem._node[2].x
                + 1.0/2 * (1-xi)/2 * elem._node[3].x;

  float dydxi = (-1.0/2) * (1-eta)/2 * elem._node[0].y
              + 1.0/2 * (1-eta)/2 * elem._node[1].y
              + 1.0/2 * (1+eta)/2 * elem._node[2].y
              + (-1.0/2) * (1+eta)/2 * elem._node[3].y;

  float dydeta = (-1.0/2) * (1-xi)/2 * elem._node[0].y
                + (-1.0/2) * (1+xi)/2 * elem._node[1].y
                + 1.0/2 * (1+xi)/2 * elem._node[2].y
                + 1.0/2 * (1-xi)/2 * elem._node[3].y;

  float jacobian = dxdxi * dydeta - dxdeta * dydxi;

  float dN1dxi = (eta-1) / 4, dN1deta = (xi-1) / 4;
  float dN2dxi = (1-eta) / 4, dN2deta = -(1+xi) / 4;
  float dN3dxi = (1+eta) / 4, dN3deta = (1+xi) / 4;
  float dN4dxi = -(1+eta) / 4, dN4deta = (1-xi) / 4;

  float dN1dx = (dN1dxi*dydeta - dN1deta*dxdeta) / jacobian; 
  float dN1dy = (-dN1dxi*dydxi + dN1deta*dxdxi) / jacobian;
  float dN2dx = (dN2dxi*dydeta - dN2deta*dxdeta) / jacobian;
  float dN2dy = (-dN2dxi*dydxi + dN2deta*dxdxi) / jacobian;
  float dN3dx = (dN3dxi*dydeta - dN3deta*dxdeta) / jacobian;
  float dN3dy = (-dN3dxi*dydxi + dN3deta*dxdxi) / jacobian;
  float dN4dx = (dN4dxi*dydeta - dN4deta*dxdeta) / jacobian;
  float dN4dy = (-dN4dxi*dydxi + dN4deta*dxdxi) / jacobian;

  for (int i=1; i<=4; ++i) {
      B[i*2-1] = 0;
      B[i*2+6] = 0;
  }
  B[0] = dN1dx; B[2] = dN2dx; B[4] = dN3dx; B[6] = dN4dx;
  B[9] = dN1dy; B[11] = dN2dy; B[13] = dN3dy; B[15] = dN4dy;
  for (int i=1; i<=4; ++i) {
      B[i*2+14] = B[i*2+7];
      B[i*2+15] = B[i*2-2];
  }

  for (int i=0; i<3; ++i) {
      for (int j=0; j<8; ++j) {
          B_t[j*3+i] = B[i*8+j];
      }
  }
}

void make_D_matrix(float *D){
  // создаю D для матрицы жестоксти
  float E = 10000000;
  float nu = 0.25;// 0.25; , 0 - чтобы были гарничные условия верны
  float tmp = E * (1-nu) / ((1+nu) * (1-2*nu));
  D[0] = tmp; D[1] = tmp * nu / (1-nu); D[2] = 0;
  D[3] = tmp * nu / (1-nu); D[4] = tmp; D[5] = 0;
  D[6] = 0; D[7] = 0; D[8] = tmp * (1 - 2*nu) / (2 * (1-nu));
}

void stiffness_matrix_local(element &elem, float *Klocal) {
  std::vector<float> quad;
  quad.push_back(-1); quad.push_back(1);
  std::vector<float> quad_w;
  quad_w.push_back(1); quad_w.push_back(1);
  float *B = new float[24];
  float *B_t = new float[24];
  float *D = new float[9];

  make_D_matrix(D);

  for (int i=0; i<64; ++i) {
      Klocal[i] = 0;
  }
  // Можно не пересчитывать много раз B, B_t
  // а сделать цикл по i, j внутри цикла по k, l.
  // Поменял циклы местами (цикл по i, j был снаружи)
  // после этого результат получился другой.
  for (int k=0; k<quad.size(); ++k) {
      for (int l=0; l<quad.size(); ++l) {
          // чтобы не пересчитывать по несолько раз
          make_grad_matrix(elem, B, B_t, quad[k], quad[l]);
          float jacobian = Jacobian(elem, quad[k], quad[l]);

          for (int i=0; i<8; ++i) {
              for (int j=0; j<8; ++j) {
                  float K_tmp = 0;

                  for (int m=0; m<3; ++m) {
                      float tmp = 0;

                      for (int n=0; n<3; ++n) {
                          tmp += D[m*3+n] * B[n*8+j];
                      }
                      K_tmp += B_t[i*3+m] * tmp;
                  }
                  K_tmp *= jacobian
                        * quad_w[k] * quad_w[l];
                  Klocal[i*8+j] += K_tmp;
              }
          }
      }
  }

  /*for (int i=0; i<8; ++i) {
	  for (int j=0; j<8; ++j) {

		for (int k=0; k<quad.size(); ++k) {
			for (int l=0; l<quad.size(); ++l) {
				make_grad_matrix(elem, B, B_t, quad[k], quad[l]);
				float jacobian = Jacobian(elem, quad[k], quad[l]);
				
				//if (jacobian < 0) 
					//std::cout << elem.eid << ": " << jacobian << std::endl;

				float K_tmp = 0;

                  for (int m=0; m<3; ++m) {
                      float tmp = 0;

                      for (int n=0; n<3; ++n) {
                          tmp += D[m*3+n] * B[n*8+j];
                      }
                      K_tmp += B_t[i*3+m] * tmp;
                  }*/
				  /*if (i==0 && j ==3) {
					  for (int r=0; r<3; ++r) {
						  for (int t=0; t<8; ++t) {
							std::cout << B[r*8+t] << "  ";
						  }
						  std::cout << std::endl;
					  }
					  std::cout << std::endl;
				  }*/
                  /*K_tmp *= jacobian
                        * quad_w[k] * quad_w[l];
				  Klocal[i*8+j] += K_tmp;
			}
		}
	  
	  }
  }*/
}

void force_matrix_local(element &elem, float *Flocal, float f, float *f_vector) {
  float ro = 10000;
  // dim of f_vector == 2
  // реализовать, но не понимаю как учитываается зависимость силы от времени
  std::vector<float> quad;
  quad.push_back(-1); quad.push_back(1);
  std::vector<float> quad_w;
  quad_w.push_back(1); quad_w.push_back(1);

  for (int i=0; i<8; ++i) {
      float tmp=0;
      // квадратуры гауса
      for (int k=0; k<quad.size(); ++k) {
          for (int l=0; l<quad.size(); ++l) {
              tmp += f * basis_function(i/2, quad[k], quad[l])
                  * Jacobian(elem, quad[k], quad[l])
                  * quad_w[k] * quad_w[l];
          }
      }
      Flocal[i] = tmp * ro * f_vector[int(i % 2)];
  }
}

#endif
