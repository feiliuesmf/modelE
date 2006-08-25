#include "rundeck_opts.h"
      SUBROUTINE DIAG_RES
!@sum  DIAG_RES Set resolution dependent diagnostic model variables 
!@sum  Coarse Resolution (36x24)
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im
      USE DIAG_COM, only : ijdd,namdd,NDIUPT
      IMPLICIT NONE

      if (im.ne.36) call stop_model
     *     ('Incorrect resolution for DIAG_RES version',255) 

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      IJDD = RESHAPE(
     &   (/31, 8,   8,17,  18,13,   6,11,  21,14,
     *     20,15,  18,15,  17,15,  19,14,  24,16,
     *     20,14,  11, 6,  17,15,  20,15,  21,15,
     *     22,16,  19,13,  18,14,  24,17,  32, 8,
     *     32, 8,  31, 8,  20,14,  18,14,  28,17,
     *     19,15,  25,15,  31, 8,  11, 7,  18,15,
     *     11, 6,  20,14,  18,14,  20,15  /),
     *     (/2,NDIUPT/) )
      NAMDD =
     &   (/'AUSD', 'MWST', 'SAHL', 'EPAC', 'AF01',
     &     'AF02', 'AF03', 'AF04', 'AF05', 'ASA1',
     &     'AF06', 'AME1', 'AF07', 'AF08', 'AF09',
     &     'ARAB', 'AF10', 'AF11', 'ASA2', 'AUS1',
     &     'AUS2', 'AUS3', 'AF12', 'AF13', 'ASA3',
     &     'AF14', 'ASA4', 'AUS4', 'AME2', 'AF15',
     &     'AME3', 'AF16', 'AF17', 'AF18' /)
#else
      IJDD = RESHAPE( (/ 31, 8,   8,17,  18,13,   6,11  /),
     *       (/2,NDIUPT/))
      NAMDD= (/ 'AUSD', 'MWST', 'SAHL', 'EPAC' /)
#endif

      RETURN
      END
