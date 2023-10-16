tryCatch(dyn.unload("src/FNKellyFunctions.dll"))
system("R CMD SHLIB src/FNKellyFunctions.f")
dyn.load("src/FNKellyFunctions.dll")

FCalcSpecConc = function(CConc, K, Stoich, NComp, NSpec){
  results = .Fortran("FNCalcSpecConc",
                     NComp = as.integer(NComp),
                     NSpec = as.integer(NSpec),
                     CConc = as.double(CConc),
                     K = as.double(K),
                     Stoich = as.integer(Stoich),
                     SConc = double(NSpec))
  return(as.numeric(results$SConc))
}
