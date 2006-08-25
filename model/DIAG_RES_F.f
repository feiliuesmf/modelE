#include "rundeck_opts.h"
      SUBROUTINE DIAG_RES
!@sum  DIAG_RES Set resolution dependent diagnostic model variables 
!@sum  Fine Resolution (144x90)
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im
      USE DIAG_COM, only : ijdd,namdd,NDIUPT
      IMPLICIT NONE

      if (im.ne.144) call stop_model
     *     ('Incorrect resolution for DIAG_RES version',255) 

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      IJDD = RESHAPE(
     &   (/125,32,  33,66,  73,52,  25,44,  85,54,
     *      81,60,  71,58,  67,58,  77,54,  97,64,
     *      79,54,  45,24,  69,58,  81,58,  83,60,
     *      89,62,  77,52,  71,56,  97,66, 127,30,
     *     127,32, 123,30,  81,54,  71,54, 113,66,
     *      75,60, 101,58, 125,30,  45,26,  73,58,
     *      45,22,  81,56,  73,56,  79,60   /),
     *       (/2,NDIUPT/))
      NAMDD =
     &   (/'AUSD', 'MWST', 'SAHL', 'EPAC', 'AF01',
     &     'AF02', 'AF03', 'AF04', 'AF05', 'ASA1',
     &     'AF06', 'AME1', 'AF07', 'AF08', 'AF09',
     &     'ARAB', 'AF10', 'AF11', 'ASA2', 'AUS1',
     &     'AUS2', 'AUS3', 'AF12', 'AF13', 'ASA3',
     &     'AF14', 'ASA4', 'AUS4', 'AME2', 'AF15',
     &     'AME3', 'AF16', 'AF17', 'AF18' /)
#else
      IJDD = RESHAPE( (/ 125,32, 33,66,  73,52,  25,44 /),
     *       (/2,NDIUPT/))
      NAMDD= (/ 'AUSD', 'MWST', 'SAHL', 'EPAC' /)
#endif

      RETURN
      END
