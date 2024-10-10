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
          Sik = SpecStoich(i, BLComp);
          dVidCj = dVdC(i, j);
          dCidCj = dCdC(i, j);
          Tmp = Sik * (dCidCj * Vi + Ci * dVidCj);
          dTkcalcdCj += Tmp;
        }
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


//' @title Calculate the Jacobian matrix
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
//' @param CompName character vector (NComp), names of components
//' @param TotConc numeric vector (NComp), input total concentrations
//' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
//'   each reaction
//' @param SpecConc numeric vector (NSpec), the concentrations of each species
//' @param SpecMC integer vector (NSpec), the index of the mass compartment to 
//'   which each species belongs
//' @param SpecCtoM numeric vector (NSpec), the factor to apply to
//'   concentrations (i.e., units of mol/L or mol/kg) to convert to masses (i.e.,
//'   units of mol).
//' @param SpecType character vector (NSpec), the type of each species
//' @param SpecCharge integer vector (NSpec), the ionic charge of the species
//' @param SpecK numeric vector (NSpec), the equilibrium coefficient of each
//'   species formation reaction
//' @param WHAMIonicStrength double, the calculated ionic strength of the 
//'   solution, excluding organic matter contribtions
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
//' @param MassAmtAdj numeric vector (NMass), the adjusted volume of each mass 
//'   compartment (adjusted for Donnan layers)
//' @param AqueousMC integer, the position in the mass compartment vectors of 
//'   the aqueous/water mass compartment
//' @param WHAMDonnanMC integer vector (2), the position in the mass compartment
//'   vectors of the DonnanHA and DonnanFA compartments
//' @param MetalComp integer, the position in component vectors of the toxic
//'   metal component
//' @param NBLMetal integer, the number of BL-Metal toxic species
//' @param BLMetalSpecs integer vector, the position in the species vectors of
//'   the metal-biotic ligand species associated with toxicity
//' @param DoTox logical, if TRUE = toxicity run, FALSE = speciation run
//'
//' @return numeric matrix (NComp x NComp), the jacobian matrix (see details)
//' @keywords internal
//'
Rcpp::NumericMatrix Jacobian1 (int NComp, //number of components
                              int NSpec, //number of species
                              Rcpp::CharacterVector CompName, //names of components
                              Rcpp::NumericVector TotConc,
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::IntegerVector SpecMC,
                              Rcpp::NumericVector SpecCtoM, //concentration to mass conversion for each species
                              Rcpp::CharacterVector SpecType,
                              Rcpp::IntegerVector SpecCharge,
                              Rcpp::NumericVector SpecK,
                              double WHAMIonicStrength,
                              bool DoWHAM,
                              Rcpp::NumericVector HumicSubstGramsPerLiter,
                              Rcpp::NumericVector WHAMSpecCharge,
                              Rcpp::NumericVector wP,
                              Rcpp::NumericVector wMolWt,
                              Rcpp::NumericVector wRadius,
                              double wDLF,
                              double wKZED,
                              Rcpp::NumericVector MassAmtAdj,
                              int AqueousMC,
                              Rcpp::IntegerVector WHAMDonnanMC,
                              int MetalComp, //position of the metal component
                              int NBLMetal, //number of BL-Metal species
                              Rcpp::IntegerVector BLMetalSpecs, //positions of BL-metal species
                              bool DoTox) {//TRUE = do toxicity, FALSE = speciation only
  
  /* outputs: */
  Rcpp::NumericMatrix JacobianMatrix(NComp);
    Rcpp::rownames(JacobianMatrix) = CompName;
    Rcpp::colnames(JacobianMatrix) = CompName;

  /* variables: */
  double Sum, Sum2;
  int iComp1, iComp2, iSpec, iHS;
  //int i;
  double Sik, Sij, Ci, Cj, Vi, Th, Hi;
  std::string Name1, Name2, NameS;
  Rcpp::String HSName;   
  Rcpp::NumericVector W = wP * log10(WHAMIonicStrength);
  Rcpp::NumericVector HumicTerm = HumicSubstGramsPerLiter * WHAMSpecCharge / Rcpp::abs(WHAMSpecCharge);
    if (WHAMSpecCharge[iHA] == 0) { HumicTerm[iHA] = 0.0; }
    if (WHAMSpecCharge[iFA] == 0) { HumicTerm[iFA] = 0.0; }
  Rcpp::NumericVector absZhThCheck = Rcpp::abs(WHAMSpecCharge) * HumicSubstGramsPerLiter;
  Rcpp::NumericMatrix dZhdC(2, NComp);
  Rcpp::NumericMatrix dVDLmaxdC(2, NComp);
  Rcpp::NumericMatrix dVDLdC(2, NComp);
  Rcpp::NumericVector MaxVolDiffusePerGramHS(2);
  Rcpp::NumericVector MaxVolDiffuse(2);
  Rcpp::NumericVector VolDiffuse(2);
  Rcpp::NumericMatrix dKidC(NSpec, NComp);
  Rcpp::NumericMatrix dVidC(NSpec, NComp);
  Rcpp::CharacterVector SpecTypeWHAM(2);
    SpecTypeWHAM[iHA] = STYPE_WHAMHA;
    SpecTypeWHAM[iFA] = STYPE_WHAMFA;

  if (DoWHAM) {
  
    CalcDonnanLayerParams(NSpec, WHAMIonicStrength, wMolWt, wRadius, wDLF, wKZED, 
                          WHAMSpecCharge, HumicSubstGramsPerLiter, 
                          MaxVolDiffusePerGramHS, MaxVolDiffuse, VolDiffuse);
  
    for (iComp2 = 0; iComp2 < NComp; iComp2++) {
      Cj = SpecConc[iComp2];
      for (iHS = 0; iHS < 2; iHS++) {
        Th = HumicSubstGramsPerLiter[iHS];
        if ((Th == 0) || (WHAMSpecCharge[iHS] == 0)) {
          dVDLmaxdC(iHS, iComp2) = 0.0;
          dZhdC(iHS, iComp2) = 0.0;
        } else {
          Sum = 0.0;
          Sum2 = 0.0;
          for (iSpec = 0; iSpec < NSpec; iSpec++) {
            Sij = SpecStoich(iSpec, iComp2);
            Hi = SpecCharge[iSpec];
            Ci = SpecConc[iSpec];
            if (SpecType[iSpec] == SpecTypeWHAM[iHS]) {
              Sum += (Sij * Ci * Hi);
              Sum2 += (Ci * pow(Hi, 2));
            }
          }
          dZhdC(iHS, iComp2) = Sum / (Th * Cj + 2 * W[iHS] * SpecConc[iComp2] * Sum2);
          dVDLmaxdC(iHS, iComp2) = (dZhdC(iHS, iComp2) * MaxVolDiffuse[iHS]) / 
                  (WHAMSpecCharge[iHS] * (1 + wKZED * std::fabs(WHAMSpecCharge[iHS])));
        }        
      }      
      // These two loops must be separate because dVDLmaxdC must be calculated for
      // both HA and FA for dVDLdC
      for (iHS = 0; iHS < 2; iHS++) {
        dVDLdC(iHS, iComp2) = 
        (dVDLmaxdC(iHS, iComp2) * wDLF * (wDLF + MaxVolDiffuse[iHA] + MaxVolDiffuse[iFA]) 
            - MaxVolDiffuse[iHS] * wDLF * (dVDLmaxdC(iHA, iComp2) + dVDLmaxdC(iFA, iComp2))) / 
          pow(wDLF + MaxVolDiffuse[iHA] + MaxVolDiffuse[iFA], 2);
      }

      for (iSpec = 0; iSpec < NSpec; iSpec++) {
        if (iSpec < NComp) {
          // components have a Ki of 1, so dKi/dCj is 0
          dKidC(iSpec, iComp2) = 0.0;
        } else if ((SpecType[iSpec] == STYPE_WHAMHA) || 
                   (SpecType[iSpec] == STYPE_WHAMFA)) {
          // WHAM species have Ki = Kint * exp(-2*w*Hi*Zh), so 
          // dKi/dCj = Ki * (-2) * w * Hi * dZh/dCj
          if (SpecType[iSpec] == STYPE_WHAMHA) { iHS = iHA; }
          if (SpecType[iSpec] == STYPE_WHAMFA) { iHS = iFA; }
          Hi = SpecCharge[iSpec];
          dKidC(iSpec, iComp2) = SpecK[iSpec] * (-2) * W[iHS] * Hi * dZhdC(iHS, iComp2);
        } else {
          // everything else (i.e., inorganic, Donnan, and BL species) have a Ki
          // that does not vary with concentration
          dKidC(iSpec, iComp2) = 0.0;
        }

        if (SpecMC[iSpec] == AqueousMC) {
          // inorganic species
          dVidC(iSpec, iComp2) = -1 * (dVDLdC(iHA, iComp2) + dVDLdC(iFA, iComp2));
        } else if ((SpecMC[iSpec] == WHAMDonnanMC[iHA]) || 
                   (SpecMC[iSpec] == WHAMDonnanMC[iFA])) {
          // Donnan Species
          if (SpecMC[iSpec] == WHAMDonnanMC[iHA]) { iHS = iHA; }
          if (SpecMC[iSpec] == WHAMDonnanMC[iFA]) { iHS = iFA; }
          dVidC(iSpec, iComp2) = dVDLdC(iHS, iComp2);
        } else {
          // WHAM and BL components
          dVidC(iSpec, iComp2) = 0.0;
        }
      }
    
    }
  } else {
    for (iComp2 = 0; iComp2 < NComp; iComp2++) {
      for (iSpec = 0; iSpec < NSpec; iSpec++) {
        dKidC(iSpec, iComp2) = 0.0;
        dVidC(iSpec, iComp2) = 0.0;
      }
    }
  }

  /* Loop through the Jacobian */
  for (iComp1 = 0; iComp1 < NComp; iComp1++) {
    Name1 = CompName(iComp1);
    for (iComp2 = 0; iComp2 < NComp; iComp2++) {
      Name2 = CompName[iComp2];
      Sum = 0;
      if (DoTox && (iComp1 == MetalComp)) {
        // Toxicity mode the metal's derivatives are relative to the CA error
        for (int i = 0; i < NBLMetal; i++) {
          iSpec = BLMetalSpecs[i];
          Sum += (SpecStoich(iSpec, iComp2) * //SpecStoich(iSpec, iComp1) *
            SpecConc[iSpec]);// * SpecCtoM[iSpec]);
        };//NEXT iSpec
        if ((SpecConc[iComp2] == 0.0) || (Sum == 0.0)) {// || (SpecCtoM[iComp2] == 0.0)
          JacobianMatrix(iComp1, iComp2) = 0.0;
        } else {
          JacobianMatrix(iComp1, iComp2) = Sum / (SpecConc[iComp2]);// * SpecCtoM(iComp2));
        }
      } else if ((SpecType[iComp1] == STYPE_WHAMHA) || 
                   (SpecType[iComp1] == STYPE_WHAMFA)) {
        // WHAM species shouldn't have their CToM adjusted, but just to be sure, 
        //   let's take it out of the calculation.
        for (iSpec = 0; iSpec < NSpec; iSpec++) {
          Sik = SpecStoich(iSpec, iComp1);
          Sij = SpecStoich(iSpec, iComp2);
          Ci = SpecConc[iSpec];
          Cj = SpecConc[iComp2];
          Sum += (Sik * Ci * (dKidC(iSpec, iComp2) / SpecK[iSpec] + (Sij / Cj)));
        };//NEXT iSpec
        JacobianMatrix(iComp1, iComp2) = Sum;
      } else if ((SpecType[iComp1] == STYPE_DONNANFA) || 
                 (SpecType[iComp1] == STYPE_DONNANHA)) {
        // Diffuse double layer residual is a function of the humic charge, 
        //   which is itself a function of component concentrations
        iHS = -1;
        if (SpecType[iComp1] == STYPE_DONNANHA) { 
          HSName = STYPE_WHAMHA; 
          iHS = iHA;
        }
        if (SpecType[iComp1] == STYPE_DONNANFA) { 
          HSName = STYPE_WHAMFA; 
          iHS = iFA;
        }
        for (iSpec = 0; iSpec < NSpec; iSpec++) {
          Sik = SpecStoich(iSpec, iComp1);
          Sij = SpecStoich(iSpec, iComp2);
          Ci = SpecConc[iSpec];
          Vi = VolDiffuse[iHS];
          Cj = SpecConc[iComp2];
          Sum += (Sik * Ci * ((dKidC(iSpec, iComp2) / SpecK[iSpec] + Sij / Cj) * Vi + dVidC(iSpec, iComp2)));
          //Sum += (Sik * Ci * (Sij * Vi / Cj + dVidC(iSpec, iComp2)));
        };//NEXT iSpec
        Sum2 = dZhdC(iHS, iComp2) * HumicTerm[iHS];
        JacobianMatrix(iComp1, iComp2) = Sum - Sum2;
      } else {
        /* All others are based on Resid = CalcTotMoles - TotMoles */
        for (iSpec = 0; iSpec < NSpec; iSpec++) {
          Sik = SpecStoich(iSpec, iComp1);
          Sij = SpecStoich(iSpec, iComp2);
          Ci = SpecConc[iSpec];
          Vi = SpecCtoM[iSpec];
          Cj = SpecConc[iComp2];
          Sum += (Sik * ((dKidC(iSpec, iComp2) * Ci / SpecK[iSpec] + Sij * Ci / Cj) * Vi + 
                  Ci * dVidC(iSpec, iComp2)));
          /*if ((Sik != 0) && (Sij != 0)) {
            CiVi = SpecConc[iSpec] * SpecCtoM[iSpec];
            Sum += (Sik * Sij * CiVi);
          }*/
        };//NEXT iSpec
        /*if ((SpecConc[iComp2] == 0.0) || (Sum == 0.0)) {// || (SpecCtoM[iComp2] == 0.0)
          JacobianMatrix(iComp1, iComp2) = 0.0;
        } else {
          JacobianMatrix(iComp1, iComp2) = Sum / (SpecConc[iComp2]);// * SpecCtoM(iComp2));
        }*/
        Sum2 = 0;//TotConc[iComp1] * dVidC(iComp1, iComp2);
        JacobianMatrix(iComp1, iComp2) = Sum - Sum2;
      }
      if (!std::isfinite(JacobianMatrix(iComp1, iComp2))) {
        /*Rcpp::Rcout << "d.Residual(" << CompName[iComp1] << 
          ") / d.Conc(" << CompName[iComp2] << ") is nan/Inf." << std::endl;*/
        
        JacobianMatrix(iComp1, iComp2) = 0.0;  
        //throw ERROR_JACOBIAN_NAN;
      }
    };//NEXT iComp2
  };//NEXT iComp1

  return JacobianMatrix;
}