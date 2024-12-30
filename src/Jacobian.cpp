// Copyright 2024 Windward Environmental LLC
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <strings.h>
#include <math.h>
#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Calculate the Jacobian matrix
//'
//' @description Calculate the Jacobian matrix of the speciation problem
//'
//' @details This function calculates the Jacobian matrix of the speciation
//'   problem. The Jacobian matrix is an (`NComp` x `NComp`) matrix of
//'   derivatives, where each row gives the derivative of the residual each total
//'   component concentration, dR[i], with respect to, in each column, the
//'   component free concentration, dX[j]. (i.e., \eqn{Z_{i,j} = \delta Resid_{i}
//'   / \delta CompConc_{j}}).
//'
//' @param NComp integer, the number of components
//' @param NSpec integer, the number of species
//' @param CompType character vector (NComp), types of components (one of the 
//'   CTYPE_ constants)
//' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
//'   each reaction
//' @param SpecConc numeric vector (NSpec), the concentrations of each species
//' @param SpecMC integer vector (NSpec), the index of the mass compartment to 
//'   which each species belongs
//' @param SpecCtoM numeric vector (NSpec), the factor to apply to
//'   concentrations (i.e., units of mol/L or mol/kg) to convert to masses (i.e.,
//'   units of mol).
//' @param SpecType character vector (NSpec), the type of each species (one of 
//'   the STYPE_ constants)
//' @param SpecCharge integer vector (NSpec), the ionic charge of the species
//' @param SpecK numeric vector (NSpec), the equilibrium coefficient of each
//'   species formation reaction
//' @param WHAMIonicStrength double, the calculated ionic strength of the 
//'   solution, excluding organic matter contributions
//' @param DoWHAM boolean, whether this simulation has WHAM components
//' @param HumicSubstGramsPerLiter numeric vector (2), the total amount of humic
//'   substances in solution, in g HS/L
//' @param WHAMSpecCharge numeric vector (2), the total charge on each type of 
//'   humic substance
//' @param wP numeric vector (2), the P parameter from WHAM, for each type of 
//'   humic substance
//' @param wMolWt numeric vector (2), the molecular weight of the humic 
//'   substances
//' @param wRadius numeric vector (2), the radius of the of humic substances
//' @param wDLF double, the double layer overlap factor from WHAM
//' @param wKZED doucle, the KZ/KZED parameter from WHAM
//' @param AqueousMC integer, the position in the mass compartment vectors of 
//'   the aqueous/water mass compartment
//' @param MetalComp integer, the position in component vectors of the toxic
//'   metal component
//' @param BLComp integer, the position in the component vectors of the biotic 
//'   ligand component associated with toxic effects
//' @param NBLMetal integer, the number of BL-Metal toxic species
//' @param BLMetalSpecs integer vector, the position in the species vectors of
//'   the metal-biotic ligand species associated with toxicity
//' @param DoTox logical, if TRUE = toxicity run, FALSE = speciation run
//'
//' @return numeric matrix (NComp x NComp), the jacobian matrix (see details)
//' @keywords internal
//'
Rcpp::NumericMatrix Jacobian(int NComp, //number of components
                              int NSpec, //number of species
                              Rcpp::CharacterVector CompType,
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::IntegerVector SpecMC,
                              Rcpp::NumericVector SpecCtoM, //concentration to mass conversion for each species
                              Rcpp::CharacterVector SpecType,
                              Rcpp::IntegerVector SpecCharge,
                              Rcpp::NumericVector SpecK,
                              Rcpp::NumericVector SpecActivityCoef,
                              double WHAMIonicStrength,
                              bool DoWHAM,
                              Rcpp::NumericVector HumicSubstGramsPerLiter,
                              Rcpp::NumericVector WHAMSpecCharge,
                              Rcpp::NumericVector wP,
                              Rcpp::NumericVector wMolWt,
                              Rcpp::NumericVector wRadius,
                              double wDLF,
                              double wKZED,
                              int AqueousMC,
                              int MetalComp, //position of the metal component
                              int BLComp,//position of BL component
                              int NBLMetal, //number of BL-Metal species
                              Rcpp::IntegerVector BLMetalSpecs, //positions of BL-metal species
                              bool DoTox,
                 bool DodVidCj,
                 bool DodVidCjDonnan,
                 bool DodKidCj,
                 bool DoGammai,
                 bool DoJacDonnan,
                 bool DoJacWHAM) {//TRUE = do toxicity, FALSE = speciation only
  
  /* outputs: */
  Rcpp::NumericMatrix dRkdCj(NComp);

  /* variables */
  int i, j, k, h, ii;
  int Sij, Sik, zi;
  double Ci, Cj, Vi, Th, Zh, Gammai, Wh;
  double Tmp;
  double dTkcalcdCj, dTkknowndCj;
  double dCidCj, dVidCj, dZhdCj;
  Rcpp::NumericMatrix dVdC(NSpec, NComp);
  Rcpp::NumericMatrix dVhdCj(2, NComp);
  Rcpp::NumericMatrix dVhmaxdCj(2, NComp);
  Rcpp::NumericMatrix dCdC(NSpec, NComp);
  double sum_CiziSij, sum_Cizi2;
  Rcpp::NumericVector w = wP * log10(WHAMIonicStrength);
  Rcpp::NumericVector ZT_absZ = HumicSubstGramsPerLiter * WHAMSpecCharge / Rcpp::abs(WHAMSpecCharge);
  Rcpp::NumericMatrix dZdC(2, NComp);
  Rcpp::CharacterVector SpecTypeWHAM(2);
  Rcpp::NumericVector VmaxHS(2);
  Rcpp::NumericVector Vmax(2);
  Rcpp::NumericVector Vh(2);
    SpecTypeWHAM[iHA] = STYPE_WHAMHA;
    SpecTypeWHAM[iFA] = STYPE_WHAMFA;

  /*bool DodVidCj = true;
  bool DodVidCjDonnan = false;
  bool DodKidCj = true;
  bool DoGammai = true;
  bool DoJacDonnan = false;
  bool DoJacWHAM = true;*/
  
  /* calculate intermediate variables */
  if (DoWHAM) {
    CalcDonnanLayerParams(NSpec, WHAMIonicStrength, wMolWt, wRadius, wDLF, wKZED, 
                          WHAMSpecCharge, HumicSubstGramsPerLiter, 
                          VmaxHS, Vmax, Vh);
  }
  for (j = 0; j < NComp; j++) {
    Cj = SpecConc(j);

    // calculate dZ(h)/dC(j)
    if (DoWHAM) {
      for (h = 0; h < 2; h++) {
        dZdC(h, j) = 0.0;
        sum_CiziSij = 0.0;
        sum_Cizi2 = 0.0;
        Th = HumicSubstGramsPerLiter(h);
        Wh = w(h);
        for (i = 0; i < NSpec; i++) {
          Ci = SpecConc(i);
          zi = SpecCharge(i);
          Sij = SpecStoich(i, j);
          if (SpecType(i) == SpecTypeWHAM(h)){
            sum_CiziSij += (Ci * zi * Sij);
            sum_Cizi2 += (Ci * pow(zi, 2));
          }        
        }//i = 0 < NSpec
        dZdC(h, j) = sum_CiziSij / (Cj * (Th + 2 * Wh * sum_Cizi2));
      }//h = 0 < 2
    } else {
      for (h = 0; h < 2; h++) {
        dZdC(h, j) = 0.0;
      }
    }

    // calculate dC(i)/dC(j)
    for (i = 0; i < NSpec; i++) {
      Ci = SpecConc(i);
      if (DoGammai) {
        Gammai = SpecActivityCoef(i);
      } else {
        Gammai = 1.0;
      }
      Sij = SpecStoich(i, j);
      if (i == j) {
        dCidCj = 1.0;
      } else if (i < NComp) {
        // i is a component, but not the same as j
        dCidCj = 0.0;
      } else if (DodKidCj && ((SpecType(i) == STYPE_WHAMHA) || (SpecType(i) == STYPE_WHAMFA))) {
        // WHAM species
        if (SpecType(i) == STYPE_WHAMHA) {
          h = iHA;
        } else if (SpecType(i) == STYPE_WHAMFA) {
          h = iFA;
        }
        zi = SpecCharge(i);
        Wh = w(h);
        dZhdCj = dZdC(h, j);
        dCidCj = (Ci / Gammai) * (-2 * Wh * zi * dZhdCj + Sij / Cj);
      } else {
        // all other species - inorganic, BL, Donnan
        dCidCj = (Ci * Sij) / (Gammai * Cj);
      }
      dCdC(i, j) = dCidCj;
    }//i = 0 < NSpec

    // calculate dV(i)/dC(j)
    if (DoWHAM && DodVidCj) {
      for (h = 0; h < 2; h++) {
        Zh = WHAMSpecCharge(h);
        dZhdCj = dZdC(h, j);
        dVhmaxdCj(h, j) = dZhdCj * Vmax(h) / (Zh * (1 + wKZED * std::abs(Zh)));
      }
      for (h = 0; h < 2; h++) {
        dVhdCj(h, j) = (
          (dVhmaxdCj(h, j) * wDLF) * (wDLF + Vmax(iHA) + Vmax(iFA)) - 
          (Vmax(h) * wDLF) * (dVhmaxdCj(iHA, j) + dVhmaxdCj(iFA, j))
        ) / pow((wDLF + Vmax(iHA) + Vmax(iFA)), 2);
      }
      for (i = 0; i < NSpec; i++) {
        if (SpecType(i) == STYPE_DONNANHA) { 
          dVdC(i, j) = dVhdCj(iHA, j);
        } else if (SpecType(i) == STYPE_DONNANFA) {
          dVdC(i, j) = dVhdCj(iFA, j);
        } else if ((SpecMC(i) == AqueousMC) && (SpecType(i) == STYPE_NORMAL)) {
          dVdC(i, j) = -1 * (dVhdCj(iHA, j) + dVhdCj(iFA, j));
        } else {
          dVdC(i, j) = 0.0;
        }
      }
    } else {
      for (i = 0; i < NSpec; i++) {
        dVdC(i, j) = 0.0;
      }
    }
  }//j = 0 < NComp

  /* calculate the jacobian */
  for (k = 0; k < NComp; k++) {
    if (CompType(k) == CTYPE_DONNANHA) {
      h = iHA;
    } else if (CompType(k) == CTYPE_DONNANFA) {
      h = iFA;
    } else {
      h = -999;
    }
    for (j = 0; j < NComp; j++) {
      Cj = SpecConc(j);
      dTkcalcdCj = 0.0;
      dTkknowndCj = 0.0;
      if (DoJacDonnan & ((CompType(k) == CTYPE_DONNANHA) || (CompType(k) == CTYPE_DONNANFA))) {
        dZhdCj = dZdC(h, j);
        dTkknowndCj = dZhdCj * ZT_absZ(h);
      }
      if ((CompType(k) == CTYPE_DONNANHA) || (CompType(k) == CTYPE_DONNANFA)) {
        if (DoJacDonnan) {
          for (i = 0; i < NSpec; i++) {
            Ci = SpecConc(i);
            Sij = SpecStoich(i, j);
            Vi = SpecCtoM(i);
            Sik = SpecStoich(i, k);
            dVidCj = dVdC(i, j);
            dCidCj = dCdC(i, j);
            if (DodVidCjDonnan) {
              Tmp = Sik * (dCidCj * Vi + Ci * dVidCj);
            } else {
              Tmp = Sik * dCidCj * Vi;
            }            
            dTkcalcdCj += Tmp;
          } 
        } else {
          dTkcalcdCj = 0.0;
        }        
      } else if ((!DoJacWHAM) && ((CompType(k) == CTYPE_WHAMHA) || (CompType(k) == CTYPE_WHAMFA))) {
        dTkcalcdCj = 0.0;
      } else if (DoTox && (k == MetalComp)) {
        for (ii = 0; ii < NBLMetal; ii++) {
          i = BLMetalSpecs(ii);
          Ci = SpecConc(i);
          Sij = SpecStoich(i, j);
          Vi = SpecCtoM(i);
          Sik = SpecStoich(i, k);
          dVidCj = dVdC(i, j);
          dCidCj = dCdC(i, j);
          if (DodVidCj) {
            Tmp = Sik * (dCidCj * Vi + Ci * dVidCj);
          } else {
            Tmp = Sik * dCidCj * Vi;
          }
          dTkcalcdCj += Tmp;
        }
        /*for (ii = 0; ii < NBLMetal; ii++) {
          i = BLMetalSpecs(ii);
          Ci = SpecConc(i);
          Sij = SpecStoich(i, j);
          Vi = SpecCtoM(i);
          Sik = SpecStoich(i, BLComp);
          dVidCj = dVdC(i, j);
          dCidCj = dCdC(i, j);
          Tmp = Sik * (dCidCj * Vi + Ci * dVidCj);
          dTkcalcdCj += Tmp;
        }*/

      } else {
        for (i = 0; i < NSpec; i++) {
          Ci = SpecConc(i);
          Sij = SpecStoich(i, j);
          Vi = SpecCtoM(i);
          Sik = SpecStoich(i, k);
          dVidCj = dVdC(i, j);
          dCidCj = dCdC(i, j);
          if (DodVidCj) {
            Tmp = Sik * (dCidCj * Vi + Ci * dVidCj);
          } else {
            Tmp = Sik * dCidCj * Vi;
          }
          dTkcalcdCj += Tmp;
        }
      }
      dRkdCj(k, j) = dTkcalcdCj - dTkknowndCj;
    }
  }
  return dRkdCj;
}
