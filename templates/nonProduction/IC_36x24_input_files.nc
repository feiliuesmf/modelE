    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m  ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_C12.D771201      ! observed init cond (atm. only)
GIC=GIC.8X10.modelE.1.nc        ! initial ground conditions
