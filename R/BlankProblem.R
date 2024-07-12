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
      ActCorr = integer(),
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
      Temp = numeric(),
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
      Test.Type = character(),
      Duration = character(),
      Lifestage = character(),
      Endpoint = character(),
      Quantifier = character(),
      References = character(),
      Miscellaneous = character()
    ),

    # WHAM parameters
    WHAM = list(
      DoWHAM = FALSE,
      WHAMVer = NA,
      WdatFile = NA,
      wDLF = as.numeric(NA),
      wKZED = as.numeric(NA),
      wP = array(as.numeric(NA), dim = 2, dimnames = list(c("HA","FA"))),
      wRadius = array(as.numeric(NA), dim = 2, dimnames = list(c("HA","FA"))),
      wMolWt = array(as.numeric(NA), dim = 2, dimnames = list(c("HA","FA")))
    ),

    Index = list(
      AqueousMCR = as.integer(NA),
      BioticLigMCR = as.integer(NA),
      WHAMDonnanMCR = array(as.integer(NA), dim = 2, dimnames = list(c("HA","FA")))
    )

  )



  return(Out)

}


#' @title Make a blank input problem list object
#'
#' @return A list object where each element in the list is an input for BLM
#'   functions.
#'
#' @keywords internal
BlankProblemList = function() {

  # assemble Output
  Out = list(

    # Counts
    NMass = 0L,
    NInLab = 0L,
    NInVar = 0L,
    NInComp = 0L,
    NDefComp = 0L,
    NComp = 0L,
    NSpec = 0L,
    NPhase = 0L,
    NSpecialDef = 0L,
    NBL = 0L,
    NMetal = 0L,
    NBLMetal = 0L,
    NCAT = 0L,

    # Mass Compartment List
    MassName = character(),
    MassAmt = numeric(),
    MassUnit = character(),
    AqueousMCR = as.integer(NA),
    BioticLigMCR = as.integer(NA),

    # Input Labels
    InLabName = character(),

    # Input Variables
    InVarName = character(),
    InVarMCName = character(),
    InVarMCR = integer(),
    InVarType = character(),

    # Input Components
    InCompName = character(),
    CompName = character(),
    CompCharge = integer(),
    CompMCName = character(),
    CompMCR = integer(),
    CompType = character(),
    CompActCorr = character(),
    CompSiteDens = numeric(),

    # Defined Components
    DefCompName = character(),
    DefCompFromNum = numeric(),
    DefCompFromVar = character(),
    DefCompCharge = integer(),
    DefCompMCName = character(),
    DefCompMCR = integer(),
    DefCompType = character(),
    DefCompActCorr = character(),
    DefCompSiteDens = numeric(),

    # Formation Reactions
    SpecName = character(),
    SpecMCName = character(),
    SpecMCR = integer(),
    SpecActCorr = character(),
    SpecNC = integer(),
    SpecCompList = matrix(data = 0L, nrow = 0, ncol = 0),
    SpecStoich = matrix(data = 0L, nrow = 0, ncol = 0),
    SpecLogK = numeric(),
    SpecDeltaH = numeric(),
    SpecTempKelvin = numeric(),
    SpecCharge = integer(),
    SpecK = numeric(),


    # Phase List
    PhaseName = character(),
    PhaseNC = integer(),
    PhaseCompList = matrix(data = 0, nrow = 0, ncol = 0),
    PhaseStoich = matrix(data = 0, nrow = 0, ncol = 0),
    PhaseLogK = numeric(),
    PhaseDeltaH = numeric(),
    PhaseTemp = numeric(),
    PhaseMoles = numeric(),

    # Special Definitions
    BLName = character(),
    BLCompR = integer(),
    MetalName = character(),
    MetalCompR = integer(),
    BLMetalName = character(),
    BLMetalSpecsR = integer(),
    DoWHAM = FALSE,

    # Critical Accumulation Table
    CATab = data.frame(Num = integer(), CA = numeric(), Species = character(),
                       Test.Type = character(), Duration = character(),
                       Lifestage = character(), Endpoint = character(),
                       Quantifier = character(), References = character(),
                       Miscellaneous = character()),

    # WHAM parameters
    WHAMDonnanMCR = array(as.integer(NA), dim =2, dimnames = list(c("HA","FA"))),
    wDLF = as.numeric(NA),
    wKZED = as.numeric(NA),
    wP = array(as.numeric(NA), dim = 2, dimnames = list(c("HA","FA"))),
    wRadius = array(as.numeric(NA), dim = 2, dimnames = list(c("HA","FA"))),
    wMolWt = array(as.numeric(NA), dim = 2, dimnames = list(c("HA","FA")))

  )



  return(Out)

}

