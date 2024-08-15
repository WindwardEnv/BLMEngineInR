#include <strings.h>
#include <math.h>
#include <Rcpp.h>
#include "CHESSFunctions.h"

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
//' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
//'   each reaction
//' @param SpecConc numeric vector (NSpec), the concentrations of each species
//' @param SpecCtoM numeric vector (NSpec), the factor to apply to
//'   concentrations (i.e., units of mol/L or mol/kg) to convert to masses (i.e.,
//'   units of mol).
//' @param SpecName character vector (NSpec), the species names
//' @param MetalComp integer, the position in component vectors of the toxic
//'   metal component
//' @param BLMetalSpecs integer vector, the position in the species vectors of
//'   the metal-biotic ligand species associated with toxicity
//' @param DoTox logical, if TRUE = toxicity run, FALSE = speciation run
//'
//' @return numeric matrix (NComp x NComp), the jacobian matrix (see details)
//' @keywords internal
//'
Rcpp::NumericMatrix Jacobian (int NComp, //number of components
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
                              double IonicStrength,
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
  Rcpp::NumericVector W = wP * log10(IonicStrength);
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

    CalcDonnanLayerParams(NSpec, IonicStrength, wMolWt, wRadius, wDLF, wKZED, 
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

      //Rcpp::Rcout << "dZhdC = [" << std::endl << dZhdC << std::endl << "]" << std::endl;
      //Rcpp::Rcout << "dVDLmaxdC = [" << std::endl << dVDLmaxdC << std::endl << "]" << std::endl;
      //Rcpp::Rcout << "dVDLdC = [" << std::endl << dVDLdC << std::endl << "]" << std::endl;

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

      //Rcpp::Rcout << "dVidC = [" << std::endl << dVidC << std::endl << "]" << std::endl;
      //Rcpp::Rcout << "dKidC = [" << std::endl << dKidC << std::endl << "]" << std::endl;
      
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
        /* WHAM species shouldn't have their CToM adjusted, but just to be sure, 
           let's take it out of the calculation. */
        if (SpecType[iComp1] == STYPE_WHAMHA) { 
          iHS = iHA;
        }
        if (SpecType[iComp1] == STYPE_WHAMFA) { 
          iHS = iFA;
        }
        for (iSpec = 0; iSpec < NSpec; iSpec++) {
          Sik = SpecStoich(iSpec, iComp1);
          Sij = SpecStoich(iSpec, iComp2);
          Ci = SpecConc[iSpec];
          Cj = SpecConc[iComp2];
          Sum += (Sik * Ci * (dKidC(iSpec, iComp2) / SpecK[iSpec] + (Sij / Cj)));
          /*if ((Sik != 0) && (Sij != 0) && SpecConc[iComp1]) {
            Sum += (Sik * Sij * SpecConc[iSpec] / SpecConc[iComp1]);
          }
          if ((Sik != 0) && (SpecCharge[iComp2]) && (dZhdC(iHS, iComp2) != 0)) {
            Sum += (-2 * W[iHS] * SpecCharge[iComp2] * 
                    Sik * SpecConc[iSpec] * dZhdC(iHS, iComp1));
          }*/
        };//NEXT iSpec
        JacobianMatrix(iComp1, iComp2) = Sum;
      } else if ((SpecType[iComp1] == STYPE_DONNANFA) || 
                 (SpecType[iComp1] == STYPE_DONNANHA)) {
        /* diffuse double layer residual is a function of the humic charge, 
           which is itself a function of component concentrations */
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
          /*if ((Sik != 0) && (Sij != 0)) {
            //term 1 = sum(S_i_DL * S_i_j * C_i / C_j) * VDLh
            Sum += (Sik * Vi * Ci * (dKidC[iSpec, iComp2] / SpecK[iSpec] + (Sij / Cj)));
          }
          if (Sik != 0) {
            //term 2 = sum(S_i_DL * C_i) * dVDLh/dCj
            Sum += (Sik * Ci * dVDLdC(iHS, iComp2));
          }*/
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