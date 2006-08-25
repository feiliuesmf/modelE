#include "rundeck_opts.h"
      SUBROUTINE DIAG_RES
!@sum  DIAG_RES Set resolution dependent diagnostic model variables 
!@sum  Medium Resolution (72x46)
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im
      USE DIAG_COM, only : ijdd,namdd,NDIUPT
      IMPLICIT NONE

      if (im.ne.72) call stop_model
     *     ('Incorrect resolution for DIAG_RES version',255) 

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      IJDD = RESHAPE(
     &   (/63,17,  17,34,  37,27,  13,23,  43,28,
     &     41,31,  36,30,  34,30,  39,28,  49,33,
     &     40,28,  23,13,  35,30,  41,30,  42,31,
     &     45,32,  39,27,  36,29,  49,34,  64,16,
     &     64,17,  62,16,  41,28,  36,28,  57,34,
     &     38,31,  51,30,  63,16,  23,14,  37,30,
     &     23,12,  41,29,  37,29,  40,31 /),
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
      IJDD = RESHAPE( (/  63,17,  17,34,  37,27,  13,23 /),
     *       (/2,NDIUPT/))
      NAMDD= (/ 'AUSD', 'MWST', 'SAHL', 'EPAC' /)
#endif

      RETURN
      END
