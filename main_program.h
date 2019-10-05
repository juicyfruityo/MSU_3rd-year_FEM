#include "load_mesh.h"
#include "make_local_matrix.h"
#include <cmath>
#include <iostream>


void main_program(float *MasMatrix, float *StifMatrix, float *ForceMatrix, int n_size) {
  load_mesh();

  // Nodes.clear();
  // делаю матрицу масс - вроде нормально работает
  float *Mloc = new float[64];
  for (int i=0; i<n_size; ++i) {
      MasMatrix[i] = 0;
  }
  std::cout << "OK" << std::endl;
  for (int i=0; i<Elements.size(); ++i) {
      mass_matrix_local(Elements[i], Mloc);
	  assembly_mass_matrix(Elements[i], MasMatrix, Mloc, n_size);

      for(int j=0; j<64; ++j) {
          Mloc[j] = 0;
      }
  }
  std::cout << "OK" << std::endl;
  for (int i=0; i<n_size; ++i) {
	  if (MasMatrix[i] == 0) {
		std::cout << "Error !!! in mass matrix:  " << i << std::endl;
		break;
	  }
  }
  std::cout << "OK" << std::endl;
  // add stiffness - need check
  float *Kloc = new float[64];
  for (int i=0; i<n_size*n_size; ++i)
      StifMatrix[i] = 0;

  for (int i=0; i<Elements.size(); ++i) {
      stiffness_matrix_local(Elements[i], Kloc);
      assembly_one_matrix(Elements[i], StifMatrix, Kloc, n_size);

      for(int j=0; j<64; ++j) {
          Kloc[j] = 0;
      }
  }
  std::cout << "OK" << std::endl;
  for (int i=0; i<n_size; ++i) {
	  if ( StifMatrix[0*n_size+i] != 0) {
		  std::cout << "StifMatrix - Ok" << std::endl;
		  break;
	  }
	  if (i == (n_size - 1)) {
		  std::cout << "Error all zeros stiff!!" << std::endl;
	  }
  }
  std::cout << "OK" << std::endl;
  //  add force - working normally
  float *Floc = new float[8];
  float *f_vector = new float[2];
  f_vector[0] = 0;  // x compoent of vector f
  f_vector[1] = -1; // y compoent of vector f
  float f = 4000; // something like f(t)
  for (int i=0; i<n_size; ++i)
      ForceMatrix[i] = 0;

  // CHANGE FORCE HERE
  // force_matrix_local(Elements[207], Floc, f);
  // assembly_force_matrix(Elements[207], ForceMatrix, Floc);
  int y_force = 10;
  /*for (int i=0; i<Elements.size(); ++i) {
	  if (Elements[i]._node[0].y == y_force || Elements[i]._node[1].y == y_force
		  || Elements[i]._node[2].y == y_force || Elements[i]._node[3].y == y_force){
			force_matrix_local(Elements[i], Floc, f);
			assembly_force_matrix(Elements[i], ForceMatrix, Floc);
			
			for (int i=0; i<8; ++i) {
				Floc[i] = 0;
			}
	  }
  }*/
  force_matrix_local(Elements[2], Floc, f, f_vector); // 50 - for lemb 10x5
  assembly_force_matrix(Elements[2], ForceMatrix, Floc);

  for (int i=0; i<n_size; ++i) {
	  if ( ForceMatrix[i] != 0) {
		  std::cout << "ForceMatrix - Ok" << std::endl;
		  break;
	  }
	  if (i == (n_size - 1)) {
		  std::cout << "Error all zeros force!!" << std::endl;
	  }
  }
  std::cout << std::endl;

  Elements.clear();

  delete[] Kloc, Floc, Mloc;

}
