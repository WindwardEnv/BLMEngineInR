#ifndef KELLYCPPSCRAP_H_INCLUDED
#define KELLYCPPSCRAP_H_INCLUDED

#include <cmath>
#include <vector>
#include <matrix.h>

vecotr<float> CppCalcSpeciesConc(unsigned int NComp, unsigned int NSpec,
                                 vector<float> CConc(NComp),
                                 vector<float> K(NSpec),
                                 Matrix<signed int> Stoich(NSpec, NComp)){
  float Tmp;
  vector<float> SConc(NSpec);

  for (unsigned int iSpec = 1; iSpec <= NSpec; iSpec++){
    Tmp = 1;
    for (unsigned int iComp = 1; iComp <= NComp; iComp++){
      Tmp = Tmp * std::pow(CConc(iComp), Stoich(iSpec, iComp));
    }
    SConc(iSpec) = Tmp
  }

  return SConc;
}


#endif // KELLYCPPSCRAP_H_INCLUDED
