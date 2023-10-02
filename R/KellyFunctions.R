CalcSpecConc = function(CConc, K, Stoich, NComp = length(CConc),
                        NSpec = length(K)){#method=2

  # if(method == 1){
  #   # Method 1
  #   SConc = rep(1, NSpec)
  #   for (iSpec in 1:NSpec){
  #     Tmp = 1
  #     for (iComp in 1:NComp){
  #       Tmp = Tmp * CConc[iComp] ^ Stoich[iSpec, iComp]
  #     }
  #     SConc[iSpec] = Tmp * K[iSpec]
  #   }
  # } else if (method == 2){
  #   # Method 2
    SConc = rep(1, NSpec)
    for (iSpec in 1:NSpec){
      X = CConc ^ Stoich[iSpec,]
      SConc[iSpec] = prod(X) * K[iSpec]
    }
  # } else if (method == 3){
  #   # Method 3
  #   SConc = apply(CConc^t(Stoich), MAR = 2, FUN = prod) * K
  # }

  return(SConc)
}

