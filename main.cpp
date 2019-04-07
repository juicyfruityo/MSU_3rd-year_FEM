#include "load_mesh.h"
#include "make_local_matrix.h"
#include <cmath>


int main() {
  load_mesh();

  Nodes.clear();

  // делаю матрицу масс - вроде нормально работает
  float *MasMatrix = new float[64*64];
  float *Mloc = new float[64];
  for (int i=0; i<64*64; ++i)
      MasMatrix[i] = 0;

  for (int i=0; i<Elements.size(); ++i) {
      mass_matrix_local(Elements[i], Mloc);
      assembly_one_matrix(Elements[i], MasMatrix, Mloc);

      for(int j=0; j<64; ++j) {
          Mloc[i] = 0;
      }
  }

  float *StifMatrix = new float[64*64];
  float *Kloc = new float[64];
  for (int i=0; i<64*64; ++i)
      StifMatrix[i] = 0;

  for (int i=0; i<Elements.size(); ++i) {
      stiffness_matrix_local(Elements[i], Kloc);
      assembly_one_matrix(Elements[i], StifMatrix, Kloc);

      for(int j=0; j<64; ++j) {
          Kloc[i] = 0;
      }
  }

  float *ForceMatrix = new float[64];
  float *Floc = new float[8];
  float f = 1; // something like f(t)
  for (int i=0; i<64; ++i)
      ForceMatrix[i] = 0;

  for (int i=0; i<Elements.size(); ++i) {
      force_matrix_local(Elements[i], Floc, f);
      assembly_force_matrix(Elements[i], ForceMatrix, Floc);

      for(int j=0; j<8; ++j) {
          Floc[i] = 0;
      }
  }


  delete[] Mloc, Kloc, Floc;
  delete[] MasMatrix, StifMatrix, ForceMatrix;



  return 0;
}
