#include "load_mesh.h"
#include "make_local_matrix.h"
#include <cmath>


int main() {
  load_mesh();

  Nodes.clear();

  // количество узлов * 2 (количество координат (x y))
  int n_size = 561 * 2;
  // делаю матрицу масс - вроде нормально работает
  float *MasMatrix = new float[n_size*n_size];
  float *Mloc = new float[n_size];
  for (int i=0; i<n_size*n_size; ++i)
      MasMatrix[i] = 0;

  for (int i=0; i<Elements.size(); ++i) {
      mass_matrix_local(Elements[i], Mloc);
      assembly_one_matrix(Elements[i], MasMatrix, Mloc, n_size);

      for(int j=0; j<n_size; ++j) {
          Mloc[i] = 0;
      }
  }

  float *StifMatrix = new float[n_size*n_size];
  float *Kloc = new float[n_size];
  for (int i=0; i<n_size*n_size; ++i)
      StifMatrix[i] = 0;

  for (int i=0; i<Elements.size(); ++i) {
      stiffness_matrix_local(Elements[i], Kloc);
      assembly_one_matrix(Elements[i], StifMatrix, Kloc, n_size);

      for(int j=0; j<n_size; ++j) {
          Kloc[i] = 0;
      }
  }

  float *ForceMatrix = new float[n_size];
  float *Floc = new float[8];
  float f = 1; // something like f(t)
  for (int i=0; i<n_size; ++i)
      ForceMatrix[i] = 0;

  // for (int i=0; i<Elements.size(); ++i) {
  force_matrix_local(Elements[273], Floc, f);
  assembly_force_matrix(Elements[273], ForceMatrix, Floc, n_size);

  //     for(int j=0; j<8; ++j) {
  //         Floc[i] = 0;
  //     }
  // }


  delete[] Mloc, Kloc, Floc;
  delete[] MasMatrix, StifMatrix, ForceMatrix;



  return 0;
}
