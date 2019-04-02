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
      mass_matrix_local(Elements[i], Kloc);
      assembly_one_matrix(Elements[i], StifMatrix, Kloc);

      for(int j=0; j<64; ++j) {
          Kloc[i] = 0;
      }
  }


  delete[] Mloc, Kloc;
  delete[] MasMatrix, StifMatrix;



  return 0;
}
