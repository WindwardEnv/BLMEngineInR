# Copyright 2024 Windward Environmental LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @title Make a blank input problem list object
#'
#' @return A list object with a template for defining the chemical problem for
#'   the `BLMEngineInR` functions. Each element in the list is a vector, list,
#'   or data.frame object grouping related parameters together. See
#'   `str(BlankProblem())` for the structure and names of the list object.
#'
#' @export
#'
#' @example tests/examples/examples-BlankProblem.R
#'
#' @family problem manipulation functions
BlankProblem = function() {

  # assemble Output
  Out = list(

    ParamFile = "",

    # Counts
    N = c(
      Mass = 0L,
      InLab = 0L,
      InVar = 0L,
      InMass = 0L,
      InComp = 0L,
      InDefComp = 0L,
      InSpec = 0L,
      DefComp = 0L,
      Comp = 0L,
      Spec = 0L,
      Phase = 0L,
      BL = 0L,
      Metal = 0L,
      BLMetal = 0L,
      CAT = 0L
    ),

    # Mass Compartment List
    Mass = data.frame(
      Name = character(),
      Amt = numeric(),
      Unit = character()
    ),


    # Input Labels
    InLabName = character(),

    # Input Variables
    InVar = data.frame(
      Name = character(),
      MCName = character(),
      MCR = integer(),
      Type = character()
    ),

    # Input components, defined components, and species
    InMassName = character(),
    InCompName = character(),
    InDefCompName = character(),
    InSpecName = character(),

    # Defined Components
    DefComp = data.frame(
      Name = character(),
      FromNum = numeric(),
      FromVar = character(),
      Charge = integer(),
      MCName = character(),
      MCR = integer(),
      Type = character(),
      ActCorr = character(),
      SiteDens = numeric()
    ),

    # Components
    Comp = data.frame(
      Name = character(),
      Charge = integer(),
      MCName = character(),
      MCR = integer(),
      Type = character(),
      ActCorr = character(),
      SiteDens = numeric()
    ),

    # Formation Reactions
    Spec = data.frame(
      Name = character(),
      Equation = character(),
      Charge = integer(),
      MCName = character(),
      MCR = integer(),
      Type = character(),
      ActCorr = character(),
      LogK = numeric(),
      K = numeric(),
      DeltaH = numeric(),
      TempKelvin = numeric(),
      NC = integer()
    ),
    SpecCompList = matrix(data = 0L, nrow = 0, ncol = 0),
    SpecStoich = matrix(data = 0L, nrow = 0, ncol = 0),

    # Phase List
    Phase = data.frame(
      Name = character(),
      Equation = character(),
      NC = integer(),
      LogK = numeric(),
      DeltaH = numeric(),
      TempKelvin = numeric(),
      Moles = numeric()
    ),
    PhaseCompList = matrix(data = 0L, nrow = 0, ncol = 0),
    PhaseStoich = matrix(data = 0L, nrow = 0, ncol = 0),


    # Special Definitions
    BL = data.frame(
      Name = character(),
      CompR = integer()
    ),
    Metal = data.frame(
      Name = character(),
      CompR = integer()
    ),
    BLMetal = data.frame(
      Name = character(),
      SpecsR = integer()
    ),

    # Critical Accumulation Table
    CATab = data.frame(
      Num = integer(),
      CA = numeric(),
      Species = character(),
      TestType = character(),
      Duration = character(),
      Lifestage = character(),
      Endpoint = character(),
      Quantifier = character(),
      References = character(),
      Miscellaneous = character()
    ),

    # WHAM parameters
    DoWHAM = FALSE,
    WHAM = BlankWHAM(),

    Index = list(
      AqueousMCR = NA_integer_,
      BioticLigMCR = NA_integer_,
      WHAMDonnanMCR = array(-1L, dim = 2, dimnames = list(c("HA", "FA")))
    )

  )

  Out

}


#' @title Make a blank input problem list object
#'
#' @return A list object where each element in the list is an input for BLM
#'   functions.
#'
#' @keywords internal
BlankProblemList = function() {

  NamedIntegerVector = integer()
  NamedNumericVector = numeric()
  NamedCharacterVector = character()
  names(NamedIntegerVector) = character()
  names(NamedNumericVector) = character()
  names(NamedCharacterVector) = character()

  WHAMNAVec = c(HA = NA_real_, FA = NA_real_)

  Out = list(

    ParamFile = "",

    # Counts
    NMass = 0L,
    NInLab = 0L,
    NInVar = 0L,
    NInMass = 0L,
    NInComp = 0L,
    NInDefComp = 0L,
    NInSpec = 0L,
    NDefComp = 0L,
    NComp = 0L,
    NSpec = 0L,
    NPhase = 0L,
    NBL = 0L,
    NMetal = 0L,
    NBLMetal = 0L,
    NCAT = 0L,

    # Mass Compartment List
    InMassName = character(),
    MassName = character(),
    MassAmt = NamedNumericVector,
    MassUnit = NamedCharacterVector,
    AqueousMCR = NA_integer_,
    BioticLigMCR = NA_integer_,

    # Input Labels
    InLabName = character(),

    # Input Variables
    InVarName = character(),
    InVarMCName = NamedCharacterVector,
    InVarMCR = NamedIntegerVector,
    InVarType = NamedCharacterVector,

    # Input Components
    InCompName = character(),
    CompName = character(),
    CompCharge = NamedIntegerVector,
    CompMCName = NamedCharacterVector,
    CompMCR = NamedIntegerVector,
    CompType = NamedCharacterVector,
    CompActCorr = NamedCharacterVector,
    CompSiteDens = NamedNumericVector,

    # Defined Components
    InDefCompName = character(),
    DefCompName = character(),
    DefCompFromNum = NamedNumericVector,
    DefCompFromVar = NamedCharacterVector,
    DefCompCharge = NamedIntegerVector,
    DefCompMCName = NamedCharacterVector,
    DefCompMCR = NamedIntegerVector,
    DefCompType = NamedCharacterVector,
    DefCompActCorr = NamedCharacterVector,
    DefCompSiteDens = NamedNumericVector,

    # Formation Reactions
    InSpecName = character(),
    SpecName = character(),
    SpecEquation = NamedCharacterVector,
    SpecCharge = NamedIntegerVector,
    SpecMCName = NamedCharacterVector,
    SpecMCR = NamedIntegerVector,
    SpecType = NamedCharacterVector,
    SpecActCorr = NamedCharacterVector,
    SpecLogK = NamedNumericVector,
    SpecK = NamedNumericVector,
    SpecDeltaH = NamedNumericVector,
    SpecTempKelvin = NamedNumericVector,
    SpecNC = NamedIntegerVector,
    SpecCompList = matrix(data = 0L, nrow = 0, ncol = 0),
    SpecStoich = matrix(data = 0L, nrow = 0, ncol = 0),

    # Phase List
    PhaseName = character(),
    PhaseEquation = NamedCharacterVector,
    PhaseNC = NamedIntegerVector,
    PhaseCompList = matrix(data = 0L, nrow = 0, ncol = 0),
    PhaseStoich = matrix(data = 0L, nrow = 0, ncol = 0),
    PhaseLogK = NamedNumericVector,
    PhaseDeltaH = NamedNumericVector,
    PhaseTempKelvin = NamedNumericVector,
    PhaseMoles = NamedNumericVector,

    # Special Definitions
    BLName = as.character(""),
    BLCompR = -1L,
    MetalName = as.character(""),
    MetalCompR = -1L,
    BLMetalName = character(),
    BLMetalSpecsR = NamedIntegerVector,

    # Critical Accumulation Table
    CATab = data.frame(Num = integer(), CA = numeric(), Species = character(),
                       TestType = character(), Duration = character(),
                       Lifestage = character(), Endpoint = character(),
                       Quantifier = character(), References = character(),
                       Miscellaneous = character()),

    # WHAM parameters
    DoWHAM = FALSE,
    WHAMDonnanMCR = array(-1L, dim = 2, dimnames = list(c("HA", "FA"))),
    WHAMVer = NA_character_,
    WHAMFile = NA_character_,
    WHAMDLF = NA_real_,
    WHAMKZED = NA_real_,
    WHAMnA = WHAMNAVec,
    WHAMpKA = WHAMNAVec,
    WHAMpKB = WHAMNAVec,
    WHAMdpKA = WHAMNAVec,
    WHAMdpKB = WHAMNAVec,
    WHAMfprB = WHAMNAVec,
    WHAMfprT = WHAMNAVec,
    WHAMdLK1A = WHAMNAVec,
    WHAMdLK1B = WHAMNAVec,
    WHAMP = WHAMNAVec,
    WHAMRadius = WHAMNAVec,
    WHAMMolWt = WHAMNAVec,
    WHAMMonodentTable = data.frame(
      S = integer(),
      AbundDenom = integer(),
      StrongWeak = character()
    ),
    WHAMBidentTable = data.frame(
      S1 = integer(),
      S2 = integer(),
      AbundDenom = integer()
    ),
    WHAMTridentTable = data.frame(
      S1 = integer(),
      S2 = integer(),
      S3 = integer(),
      AbundDenom = integer()
    ),
    WHAMMetalsTable = data.frame(
      Metal = character(),
      pKMAHA = numeric(),
      pKMAFA = numeric(),
      dLK2 = numeric()
    ),
    WHAMSpecKselTable = data.frame(
      Spec = character(),
      KselHA = numeric(),
      KselFA = numeric()
    )
  )

  Out

}
