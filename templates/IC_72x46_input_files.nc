    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m  ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_M20A.D771201          ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext.nc ! initial conditions (ground)
