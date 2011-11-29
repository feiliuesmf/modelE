module MolecularWeights_mod
  use Dictionary_mod
  implicit none
  private

  public :: getMolecularWeight
  public :: initializeMolecularWeights

  type (Dictionary_type) :: molecularWeights ! note - private table

contains

  subroutine initializeMolecularWeights()
    use Constant, only: mair, mwat

    call add('Air', mair)
    call add('CO2n', 44.d0)
    call add('CFCn',  137.37d0)
    call add('SF6', 146.01d0)
    call add('Rn222', 222.d0)
    call add('CO2', 44.d0)
    call add('N2O', 44.d0)
    call add('CFC11',  137.4d0)
    call add('14CO2',  46.d0)
    call add('CH4', 16.d0)
    call add('O3', 48.d0)
    call add('SF6_c', 146.01d0)
    call add('Water', mwat)
    call add('H2O18', 20d0)
    call add('HDO', 19d0)
    call add('HTO', 20d0)
    call add('H2O17', 19d0)
    call add('Ox', 48.d0)
    call add('NOx', 14.01d0)
    call add('N2O5', 108.02d0)
    call add('ClOx', 51.5d0)
    call add('BrOx', 95.9d0)
    call add('HCl', 36.5d0)
    call add('ClONO2', 97.5d0)
    call add('HOCl', 52.5d0)
    call add('HBr', 80.9d0)
    call add('HOBr', 96.9d0)
    call add('BrONO2', 141.9d0)
    call add('CFC', 137.4d0) !CFC11
    call add('HNO3', 63.018d0)
    call add('H2O2', 34.016d0)
    call add('GLT', mair) ! generic linear tracer
    call add('stratOx', 48.d0)
    call add('codirect', 28.01d0)
    call add('CH3OOH', 48.042d0)
    call add('HCHO', 30.026d0)
    call add('HO2NO2', 79.018d0)
    call add('CO', 28.01d0)
    call add('PAN', 121.054d0)   ! assuming CH3COOONO2 = PAN)
    call add('Isoprene', 60.05d0) ! i.e. 5 carbons
    call add('AlkylNit', mair) !unknown molecular weight,  so use air and make
    call add('Alkenes', 1.0d0)
    call add('Paraffin', 1.0d0)
    call add('Terpenes', 120.10d0) ! i.e. 10 carbons
    call add('isopp1g', 15.6d0)
    call add('isopp1a', 15.6d0)
    call add('isopp2g', 15.6d0)
    call add('isopp2a', 15.6d0)
    call add('apinp1g', 15.6d0)
    call add('apinp1a', 15.6d0)
    call add('apinp2g', 15.6d0)
    call add('apinp2a', 15.6d0)
    call add('DMS', 62.d+0)
    call add('MSA', 96.d+0)
    call add('SO2', 64.d+0)
    call add('SO4', 96.d+0)
    call add('SO4_d1', 96.d0)
    call add('SO4_d2', 96.d0)
    call add('SO4_d3', 96.d0)
    call add('N_d1', 62.d+0)
    call add('N_d2', 62.d+0)
    call add('N_d3', 62.d0)
    call add('BCII', 12.d0)  !Insoluble industrial BC
    call add('BCIA', 12.d0)  !Soluble (aged) industrial BC
    call add('BCB', 12.d0)   !Biomass BC
    call add('OCII', 15.6d0) !Insoluble industrial organic mass
    call add('OCIA', 15.6d0) !Aged industrial organic mass
    call add('OCB', 15.6d0)  !Biomass organic mass

  contains

    subroutine add(name, weight)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: weight

      call insert(MolecularWeights, name, weight)

    end subroutine add

  end subroutine initializeMolecularWeights
  
  real*8 function getMolecularWeight(name) result(weight)
    character(len=*), intent(in) :: name
    weight = lookup(molecularWeights, name)
  end function getMolecularWeight

end module MolecularWeights_mod


