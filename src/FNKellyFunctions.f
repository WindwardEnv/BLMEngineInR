      subroutine FNCalcSpecConc(NComp, NSpec, CConc, K, Stoich, SConc)

            use, intrinsic :: iso_fortran_env, only: dp=>real64

            integer, intent(in) :: NComp, NSpec
            real(dp), intent(in) :: CConc(NComp)
            real(dp), intent(in) :: K(NSpec)
            integer, intent(in) :: Stoich(NSpec, NComp)
            real(dp), intent(out) :: SConc(NSpec)
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

      end subroutine FNCalcSpecConc
