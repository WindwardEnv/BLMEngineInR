test_that("ConvertWHAM7Files works", {

  DB7File = file.path("C:","Users","Public","Documents","WHAM7 Databases","Backups","default.db7")
  HAph7File = file.path("C:","Users","Public","Documents","WHAM7 Databases","Backups","humic acid.ph7")
  FAph7File = file.path("C:","Users","Public","Documents","WHAM7 Databases","Backups","fulvic acid.ph7")

  skip_if_not(file.exists(DB7File))

  tf_param = withr::local_tempfile()
  tf_wham = withr::local_tempfile()

  expect_no_error({
    thisProb1 = ConvertWHAM7Files(
      DB7File = DB7File,
      HAph7File = HAph7File,
      FAph7File = FAph7File,
      RParamFile = tf_param,
      RWHAMFile = tf_wham
    )
  })

  expect_no_error({thisWHAM = DefineWHAM(WHAMFile = tf_wham)})
  expect_no_error({thisProb = DefineProblem(ParamFile = tf_param)})

  CompareNames = setdiff(names(thisWHAM), c("Ver", "File"))
  expect_equal(thisWHAM[CompareNames], thisProb1$WHAM[CompareNames])
  expect_equal(thisWHAM[CompareNames], thisProb$WHAM[CompareNames])

  CompareNames = setdiff(names(thisProb), c("ParamFile", "WHAM"))
  expect_equal(thisProb1[CompareNames], thisProb[CompareNames])

  # thisWHAM is WHAM 7
  expect_equal(thisWHAM$MonodentTable$S, 1:8)
  expect_equal(thisWHAM$MonodentTable$AbundDenom, rep(c(4,8), each = 4))
  expect_equal(thisWHAM$MonodentTable$StrongWeak, rep(c("S","W"), each = 4))
  expect_setequal(paste0(thisWHAM$BidentTable$S1, thisWHAM$BidentTable$S2),
                  paste0(c(1,3,1,2,3,4), c(2,4,5,6,7,8)))
  expect_equal(thisWHAM$BidentTable$AbundDenom, rep(8, each = 6))
  expect_setequal(paste0(thisWHAM$TridentTable$S1, thisWHAM$TridentTable$S2, thisWHAM$TridentTable$S3),
                  paste0(rep(c(1,3), each = 4), rep(c(2,4), each = 4), rep(c(5,6,7,8), times = 2)))
  expect_equal(thisWHAM$TridentTable$AbundDenom, rep(16, each = 8))

})
