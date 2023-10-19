      function FNCalcSpecConc(NComp, NSpec, CConc, K, Stoich) result(SConc)

          use, intrinsic :: iso_c_binding, only: dp=>c_double

          integer, intent(in) :: NComp, NSpec
          real(dp), intent(in) :: CConc(NComp)
          real(dp), intent(in) :: K(NSpec)
          integer, intent(in) :: Stoich(NSpec, NComp)
          real(dp) :: SConc(NSpec)
          integer :: iSpec
          integer :: iComp
          real(dp) :: Tmp

          do iSpec = 1, NSpec
              Tmp = 1
              do iComp = 1, NComp
                  Tmp = Tmp * (CConc(iComp) ** Stoich(iSpec, iComp))
              end do
              SConc(iSpec) = Tmp * K(iSpec)
          end do

      end function FNCalcSpecConc
