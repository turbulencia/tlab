
! Do loop at least 1 time
  CHW_ERROR = C_1EM3_R

  CHW_INC = 2 * CHW_ERROR

  DO WHILE ( ABS(CHW_INC) .GT. CHW_ERROR ) 

#include "dns_chem_enth.h"
  CHW_INC      = (MACRO_HINPUT - CH_H)/CHW_CPTOTAL
  MACRO_TINPUT = CHW_INC + MACRO_TINPUT
  
  ENDDO

