      subroutine qng(f,p,a,b,epsabs,epsrel,result,abserr,neval,ier)

      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_50_L       50.0e+0
#define C_05EM14_L   0.5e-14
#define C_200_L      200.0e+0
#define C_1P5_L      1.5e+0


#define C_X1C1_R  0.9739065285171717e+00
#define C_X1C2_R  0.8650633666889845e+00
#define C_X1C3_R  0.6794095682990244e+00
#define C_X1C4_R  0.4333953941292472e+00
#define C_X1C5_R  0.1488743389816312e+00

#define C_X2C1_R  0.9956571630258081e+00
#define C_X2C2_R  0.9301574913557082e+00
#define C_X2C3_R  0.7808177265864169e+00
#define C_X2C4_R  0.5627571346686047e+00
#define C_X2C5_R  0.2943928627014602e+00

#define C_X3C1_R  0.9993333609019321e+00
#define C_X3C2_R  0.9874334029080889e+00
#define C_X3C3_R  0.9548079348142663e+00
#define C_X3C4_R  0.9001486957483283e+00
#define C_X3C5_R  0.8251983149831142e+00
#define C_X3C6_R  0.7321483889893050e+00
#define C_X3C7_R  0.6228479705377252e+00
#define C_X3C8_R  0.4994795740710565e+00
#define C_X3C9_R  0.3649016613465808e+00
#define C_X3C10_R 0.2222549197766013e+00
#define C_X3C11_R 0.7465061746138332e-01


#define C_X4C1_R   0.9999029772627292e+00
#define C_X4C2_R   0.9979898959866787e+00 
#define C_X4C3_R   0.9921754978606872e+00
#define C_X4C4_R   0.9813581635727128e+00  
#define C_X4C5_R   0.9650576238583846e+00
#define C_X4C6_R   0.9431676131336706e+00  
#define C_X4C7_R   0.9158064146855072e+00
#define C_X4C8_R   0.8832216577713165e+00  
#define C_X4C9_R   0.8457107484624157e+00
#define C_X4C10_R  0.8035576580352310e+00  
#define C_X4C11_R  0.7570057306854956e+00
#define C_X4C12_R  0.7062732097873218e+00  
#define C_X4C13_R  0.6515894665011779e+00
#define C_X4C14_R  0.5932233740579611e+00  
#define C_X4C15_R  0.5314936059708319e+00
#define C_X4C16_R  0.4667636230420228e+00 
#define C_X4C17_R  0.3994248478592188e+00
#define C_X4C18_R  0.3298748771061883e+00 
#define C_X4C19_R  0.2585035592021616e+00
#define C_X4C20_R  0.1856953965683467e+00 
#define C_X4C21_R  0.1118422131799075e+00
#define C_X4C22_R  0.3735212339461987e-01


#define C_W10C1    0.6667134430868814e-01
#define C_W10C2    0.1494513491505806e+00
#define C_W10C3    0.2190863625159820e+00
#define C_W10C4    0.2692667193099964e+00
#define C_W10C5    0.2955242247147529e+00

#define C_W21AC1   0.3255816230796473e-01
#define C_W21AC2   0.7503967481091995e-01
#define C_W21AC3   0.1093871588022976e+00
#define C_W21AC4   0.1347092173114733e+00
#define C_W21AC5   0.1477391049013385e+00

#define C_W21BC1   0.1169463886737187e-01
#define C_W21BC2   0.5475589657435200e-01
#define C_W21BC3   0.9312545458369761e-01
#define C_W21BC4   0.1234919762620659e+00
#define C_W21BC5   0.1427759385770601e+00
#define C_W21BC6   0.1494455540029169e+00

#define C_W43AC1   0.1629673428966656e-01
#define C_W43AC2   0.3752287612086950e-01
#define C_W43AC3   0.5469490205825544e-01
#define C_W43AC4   0.6735541460947809e-01
#define C_W43AC5   0.7387019963239395e-01
#define C_W43AC6   0.5768556059769796e-02
#define C_W43AC7   0.2737189059324884e-01
#define C_W43AC8   0.4656082691042883e-01
#define C_W43AC9   0.6174499520144256e-01
#define C_W43AC10  0.7138726726869340e-01


#define C_W43BC1     0.1844477640212414e-02
#define C_W43BC2     0.1079868958589165e-01
#define C_W43BC3     0.2189536386779543e-01
#define C_W43BC4     0.3259746397534569e-01
#define C_W43BC5     0.4216313793519181e-01
#define C_W43BC6     0.5074193960018458e-01
#define C_W43BC7     0.5837939554261925e-01
#define C_W43BC8     0.6474640495144589e-01
#define C_W43BC9     0.6956619791235648e-01
#define C_W43BC10    0.7282444147183321e-01
#define C_W43BC11    0.7450775101417512e-01
#define C_W43BC12    0.7472214751740301e-01

     
#define C_W87AC1     0.8148377384149173e-02
#define C_W87AC2     0.1876143820156282e-01
#define C_W87AC3     0.2734745105005229e-01
#define C_W87AC4     0.3367770731163793e-01
#define C_W87AC5     0.3693509982042791e-01
#define C_W87AC6     0.2884872430211531e-02
#define C_W87AC7     0.1368594602271270e-01
#define C_W87AC8     0.2328041350288831e-01
#define C_W87AC9     0.3087249761171336e-01
#define C_W87AC10    0.3569363363941877e-01
#define C_W87AC11    0.9152833452022414e-03
#define C_W87AC12    0.5399280219300471e-02
#define C_W87AC13    0.1094767960111893e-01
#define C_W87AC14    0.1629873169678734e-01
#define C_W87AC15    0.2108156888920384e-01
#define C_W87AC16    0.2537096976925383e-01
#define C_W87AC17    0.2918969775647575e-01
#define C_W87AC18    0.3237320246720279e-01
#define C_W87AC19    0.3478309895036514e-01 
#define C_W87AC20    0.3641222073135179e-01
#define C_W87AC21    0.3725387550304771e-01


#define C_W87BC1     0.2741455637620724e-03
#define C_W87BC2     0.1807124155057943e-02 
#define C_W87BC3     0.4096869282759165e-02
#define C_W87BC4     0.6758290051847379e-02
#define C_W87BC5     0.9549957672201647e-02
#define C_W87BC6     0.1232944765224485e-01
#define C_W87BC7     0.1501044734638895e-01
#define C_W87BC8     0.1754896798624319e-01
#define C_W87BC9     0.1993803778644089e-01
#define C_W87BC10    0.2219493596101229e-01
#define C_W87BC11    0.2433914712600081e-01
#define C_W87BC12    0.2637450541483921e-01
#define C_W87BC13    0.2828691078877120e-01
#define C_W87BC14    0.3005258112809270e-01
#define C_W87BC15    0.3164675137143993e-01
#define C_W87BC16    0.3305041341997850e-01
#define C_W87BC17    0.3425509970422606e-01
#define C_W87BC18    0.3526241266015668e-01
#define C_W87BC19    0.3607698962288870e-01
#define C_W87BC20    0.3669860449845609e-01
#define C_W87BC21    0.3712054926983258e-01
#define C_W87BC22    0.3733422875193504e-01
#define C_W87BC23    0.3736107376267902e-01

#else
#define C_50_L       50.0d+0
#define C_05EM14_L   0.5d-14
#define C_200_L      200.0d+0
#define C_1P5_L      1.5d+0

#define C_X1C1_R  0.9739065285171717d+00
#define C_X1C2_R  0.8650633666889845d+00
#define C_X1C3_R  0.6794095682990244d+00
#define C_X1C4_R  0.4333953941292472d+00
#define C_X1C5_R  0.1488743389816312d+00

#define C_X2C1_R  0.9956571630258081d+00
#define C_X2C2_R  0.9301574913557082d+00
#define C_X2C3_R  0.7808177265864169d+00
#define C_X2C4_R  0.5627571346686047d+00
#define C_X2C5_R  0.2943928627014602d+00

#define C_X3C1_R  0.9993333609019321d+00
#define C_X3C2_R  0.9874334029080889d+00
#define C_X3C3_R  0.9548079348142663d+00
#define C_X3C4_R  0.9001486957483283d+00
#define C_X3C5_R  0.8251983149831142d+00
#define C_X3C6_R  0.7321483889893050d+00
#define C_X3C7_R  0.6228479705377252d+00
#define C_X3C8_R  0.4994795740710565d+00
#define C_X3C9_R  0.3649016613465808d+00
#define C_X3C10_R 0.2222549197766013d+00
#define C_X3C11_R 0.7465061746138332d-01


#define C_X4C1_R   0.9999029772627292d+00
#define C_X4C2_R   0.9979898959866787d+00 
#define C_X4C3_R   0.9921754978606872d+00
#define C_X4C4_R   0.9813581635727128d+00  
#define C_X4C5_R   0.9650576238583846d+00
#define C_X4C6_R   0.9431676131336706d+00  
#define C_X4C7_R   0.9158064146855072d+00
#define C_X4C8_R   0.8832216577713165d+00  
#define C_X4C9_R   0.8457107484624157d+00
#define C_X4C10_R  0.8035576580352310d+00  
#define C_X4C11_R  0.7570057306854956d+00
#define C_X4C12_R  0.7062732097873218d+00  
#define C_X4C13_R  0.6515894665011779d+00
#define C_X4C14_R  0.5932233740579611d+00  
#define C_X4C15_R  0.5314936059708319d+00
#define C_X4C16_R  0.4667636230420228d+00 
#define C_X4C17_R  0.3994248478592188d+00
#define C_X4C18_R  0.3298748771061883d+00 
#define C_X4C19_R  0.2585035592021616d+00
#define C_X4C20_R  0.1856953965683467d+00 
#define C_X4C21_R  0.1118422131799075d+00
#define C_X4C22_R  0.3735212339461987d-01


#define C_W10C1    0.6667134430868814d-01
#define C_W10C2    0.1494513491505806d+00
#define C_W10C3    0.2190863625159820d+00
#define C_W10C4    0.2692667193099964d+00
#define C_W10C5    0.2955242247147529d+00

#define C_W21AC1   0.3255816230796473d-01
#define C_W21AC2   0.7503967481091995d-01
#define C_W21AC3   0.1093871588022976d+00
#define C_W21AC4   0.1347092173114733d+00
#define C_W21AC5   0.1477391049013385d+00

#define C_W21BC1   0.1169463886737187d-01
#define C_W21BC2   0.5475589657435200d-01
#define C_W21BC3   0.9312545458369761d-01
#define C_W21BC4   0.1234919762620659d+00
#define C_W21BC5   0.1427759385770601d+00
#define C_W21BC6   0.1494455540029169d+00

#define C_W43AC1   0.1629673428966656d-01
#define C_W43AC2   0.3752287612086950d-01
#define C_W43AC3   0.5469490205825544d-01
#define C_W43AC4   0.6735541460947809d-01
#define C_W43AC5   0.7387019963239395d-01
#define C_W43AC6   0.5768556059769796d-02
#define C_W43AC7   0.2737189059324884d-01
#define C_W43AC8   0.4656082691042883d-01
#define C_W43AC9   0.6174499520144256d-01
#define C_W43AC10  0.7138726726869340d-01


#define C_W43BC1     0.1844477640212414d-02
#define C_W43BC2     0.1079868958589165d-01
#define C_W43BC3     0.2189536386779543d-01
#define C_W43BC4     0.3259746397534569d-01
#define C_W43BC5     0.4216313793519181d-01
#define C_W43BC6     0.5074193960018458d-01
#define C_W43BC7     0.5837939554261925d-01
#define C_W43BC8     0.6474640495144589d-01
#define C_W43BC9     0.6956619791235648d-01
#define C_W43BC10    0.7282444147183321d-01
#define C_W43BC11    0.7450775101417512d-01
#define C_W43BC12    0.7472214751740301d-01

     
#define C_W87AC1     0.8148377384149173d-02
#define C_W87AC2     0.1876143820156282d-01
#define C_W87AC3     0.2734745105005229d-01
#define C_W87AC4     0.3367770731163793d-01
#define C_W87AC5     0.3693509982042791d-01
#define C_W87AC6     0.2884872430211531d-02
#define C_W87AC7     0.1368594602271270d-01
#define C_W87AC8     0.2328041350288831d-01
#define C_W87AC9     0.3087249761171336d-01
#define C_W87AC10    0.3569363363941877d-01
#define C_W87AC11    0.9152833452022414d-03
#define C_W87AC12    0.5399280219300471d-02
#define C_W87AC13    0.1094767960111893d-01
#define C_W87AC14    0.1629873169678734d-01
#define C_W87AC15    0.2108156888920384d-01
#define C_W87AC16    0.2537096976925383d-01
#define C_W87AC17    0.2918969775647575d-01
#define C_W87AC18    0.3237320246720279d-01
#define C_W87AC19    0.3478309895036514d-01 
#define C_W87AC20    0.3641222073135179d-01
#define C_W87AC21    0.3725387550304771d-01


#define C_W87BC1     0.2741455637620724d-03
#define C_W87BC2     0.1807124155057943d-02 
#define C_W87BC3     0.4096869282759165d-02
#define C_W87BC4     0.6758290051847379d-02
#define C_W87BC5     0.9549957672201647d-02
#define C_W87BC6     0.1232944765224485d-01
#define C_W87BC7     0.1501044734638895d-01
#define C_W87BC8     0.1754896798624319d-01
#define C_W87BC9     0.1993803778644089d-01
#define C_W87BC10    0.2219493596101229d-01
#define C_W87BC11    0.2433914712600081d-01
#define C_W87BC12    0.2637450541483921d-01
#define C_W87BC13    0.2828691078877120d-01
#define C_W87BC14    0.3005258112809270d-01
#define C_W87BC15    0.3164675137143993d-01
#define C_W87BC16    0.3305041341997850d-01
#define C_W87BC17    0.3425509970422606d-01
#define C_W87BC18    0.3526241266015668d-01
#define C_W87BC19    0.3607698962288870d-01
#define C_W87BC20    0.3669860449845609d-01
#define C_W87BC21    0.3712054926983258d-01
#define C_W87BC22    0.3733422875193504d-01
#define C_W87BC23    0.3736107376267902d-01

#endif

!**begin prologue  qng
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
!***keywords  automatic integrator, smooth integrand,
!     non-adaptive, gauss-kronrod(patterson)
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!   de doncker,elise,appl math & progr. div. - k.u.leuven
!   kahaner,david,nbs - modified (2/82)
!***purpose  the routine calculates an approximation result to a
!    given definite integral i = integral of f over (a,b),
!    hopefully satisfying following claim for accuracy
!    abs(i-result).le.max(epsabs,epsrel*abs(i)).
!***description
!
! non-adaptive integration
! standard fortran subroutine
! real version
!
!   f      - real version
!            function subprogram defining the integrand function
!            f(x). the actual name for f needs to be declared
!            e x t e r n a l in the driver program.
!
!   a      - real version
!            lower limit of integration
!
!   b      - real version
!            upper limit of integration
!
!   epsabs - real
!            absolute accuracy requested
!   epsrel - real
!            relative accuracy requested
!            if  epsabs.le.0
!            and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!            the routine will end with ier = 6.
!
! on return
!   result - real
!            approximation to the integral i
!            result is obtained by applying the 21-point
!            gauss-kronrod rule (res21) obtained by optimal
!            addition of abscissae to the 10-point gauss rule
!            (res10), or by applying the 43-point rule (res43)
!            obtained by optimal addition of abscissae to the
!            21-point gauss-kronrod rule, or by applying the
!            87-point rule (res87) obtained by optimal addition
!            of abscissae to the 43-point rule.
!
!   abserr - real
!            estimate of the modulus of the absolute error,
!            which should equal or exceed abs(i-result)
!
!   neval  - integer
!            number of integrand evaluations
!
!   ier    - ier = 0 normal and reliable termination of the
!                    routine. it is assumed that the requested
!                    accuracy has been achieved.
!            ier.gt.0 abnormal termination of the routine. it is
!                    assumed that the requested accuracy has
!                    not been achieved.
!   error messages
!            ier = 1 the maximum number of steps has been
!                    executed. the integral is probably too
!                    difficult to be calculated by dqng.
!                = 6 the input is invalid, because
!                    epsabs.le.0 and
!                    epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
!                    result, abserr and neval are set to zero.
!
!***references  (none)
!***routines called  d1mach,xerror
!***end prologue  qng
!
      TREAL a,absc,abserr,b,centr,dhlgth,epmach,epsabs,epsrel,f,&
           fcentr,fval,fval1,fval2,fv1,fv2,fv3,fv4,hlgth,result,res10,&
           res21,res43,res87,resabs,resasc,reskh,M_MACH,&
           savfun,uflow,&
           w10,w21a,w21b, w43a,w43b,w87a,w87b,x1,x2,x3,x4,p
      TINTEGER ier,ipx,k,l,neval
      TINTEGER i1,i4
      external f
!
      dimension fv1(5),fv2(5),fv3(5),fv4(5),x1(5),x2(5),x3(11),x4(22),&
        w10(5),w21a(5),w21b(6),w43a(10),w43b(12),w87a(21),w87b(23),&
        savfun(21)
!
!   the following data statements contain the
!   abscissae and weights of the integration rules used.
!
!   x1      abscissae common to the 10-, 21-, 43-
!           and 87-point rule
!   x2      abscissae common to the 21-, 43- and
!           87-point rule
!   x3      abscissae common to the 43- and 87-point
!           rule
!   x4      abscissae of the 87-point rule
!   w10     weights of the 10-point formula
!   w21a    weights of the 21-point formula for
!           abscissae x1
!   w21b    weights of the 21-point formula for
!           abscissae x2
!   w43a    weights of the 43-point formula for
!           abscissae x1, x3
!   w43b    weights of the 43-point formula for
!           abscissae x3
!   w87a    weights of the 87-point formula for
!           abscissae x1, x2, x3
!   w87b    weights of the 87-point formula for
!           abscissae x4
!
      data x1(1),x1(2),x1(3),x1(4),x1(5)/&
           C_X1C1_R,     C_X1C2_R,&
           C_X1C3_R,     C_X1C4_R,&
           C_X1C5_R/
      data x2(1),x2(2),x2(3),x2(4),x2(5)/&
           C_X2C1_R,     C_X2C2_R,&
           C_X2C3_R,     C_X2C4_R,&
           C_X2C5_R/
      data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),&
        x3(9),x3(10),x3(11)/&
           C_X3C1_R,     C_X3C2_R,&
           C_X3C3_R,     C_X3C4_R,&
           C_X3C5_R,     C_X3C6_R,&
           C_X3C7_R,     C_X3C8_R,&
           C_X3C9_R,     C_X3C10_R,&
           C_X3C11_R/
      data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),&
        x4(10),x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),&
        x4(19),x4(20),x4(21),x4(22)/   C_X4C1_R,&
           C_X4C2_R,     C_X4C3_R,&
           C_X4C4_R,     C_X4C5_R,&
           C_X4C6_R,     C_X4C7_R,&
           C_X4C8_R,     C_X4C9_R,&
           C_X4C10_R,     C_X4C11_R,&
           C_X4C12_R,     C_X4C13_R,&
           C_X4C14_R,     C_X4C15_R,&
           C_X4C16_R,     C_X4C17_R,&
           C_X4C18_R,     C_X4C19_R,&
           C_X4C20_R,     C_X4C21_R,&
           C_X4C22_R/
      data w10(1),w10(2),w10(3),w10(4),w10(5)/&
           C_W10C1,     C_W10C2,&
           C_W10C3,     C_W10C4,&
           C_W10C5/
      data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/&
           C_W21AC1,     C_W21AC2,&
           C_W21AC3,     C_W21AC4,&
           C_W21AC5/
      data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/&
           C_W21BC1,     C_W21BC2,&
           C_W21BC3,     C_W21BC4,&
           C_W21BC5,     C_W21BC6/
      data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7),&
        w43a(8),w43a(9),w43a(10)/      C_W43AC1,&
           C_W43AC2,     C_W43AC3,&
           C_W43AC4,     C_W43AC5,&
           C_W43AC6,     C_W43AC7,&
           C_W43AC8,     C_W43AC9,&
           C_W43AC10/
      data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),&
        w43b(7),w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/&
           C_W43BC1,     C_W43BC2,&
           C_W43BC3,     C_W43BC4,&
           C_W43BC5,     C_W43BC6,&
           C_W43BC7,     C_W43BC8,&
           C_W43BC9,     C_W43BC10,&
           C_W43BC11,     C_W43BC12/
      data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),&
        w87a(7),w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),&
        w87a(13),w87a(14),w87a(15),w87a(16),w87a(17),w87a(18),&
        w87a(19),w87a(20),w87a(21)/&
           C_W87AC1,     C_W87AC2,&
           C_W87AC3,     C_W87AC4,&
           C_W87AC5,     C_W87AC6,&
           C_W87AC7,     C_W87AC8,&
           C_W87AC9,     C_W87AC10,&
           C_W87AC11,     C_W87AC12,&
           C_W87AC13,     C_W87AC14,&
           C_W87AC15,     C_W87AC16,&
           C_W87AC17,     C_W87AC18,&
           C_W87AC19,     C_W87AC20,&
           C_W87AC21/
      data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7),&
        w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14),&
        w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),&
        w87b(21),w87b(22),w87b(23)/    C_W87BC1,&
           C_W87BC2,     C_W87BC3,&
           C_W87BC4,     C_W87BC5,&
           C_W87BC6,     C_W87BC7,&
           C_W87BC8,     C_W87BC9,&
           C_W87BC10,     C_W87BC11,&
           C_W87BC12,     C_W87BC13,&
           C_W87BC14,     C_W87BC15,&
           C_W87BC16,     C_W87BC17,&
           C_W87BC18,     C_W87BC19,&
           C_W87BC20,     C_W87BC21,&
           C_W87BC22,     C_W87BC23/
!
!   list of major variables
!   -----------------------
!
!   centr  - mid point of the integration interval
!   hlgth  - half-length of the integration interval
!   fcentr - function value at mid point
!   absc   - abscissa
!   fval   - function value
!   savfun - array of function values which
!            have already been computed
!   res10  - 10-point gauss result
!   res21  - 21-point kronrod result
!   res43  - 43-point result
!   res87  - 87-point result
!   resabs - approximation to the integral of abs(f)
!   resasc - approximation to the integral of abs(f-i/(b-a))
!
!   machine dependent constants
!   ---------------------------
!
!   epmach is the largest relative spacing.
!   uflow is the smallest positive magnitude.
!
!***first executable statement  qng
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
!   test on validity of parameters
!   ------------------------------
!
      result = C_0_R
      abserr = C_0_R
      neval = 0
      ier = 6
      if(epsabs.le.C_0_R.and.epsrel.lt.&
           MAX(C_05EM14_L,C_50_L*epmach)) go to 80
      hlgth = C_05_R*(b-a)
      dhlgth = abs(hlgth)
      centr = C_05_R*(b+a)
      fcentr = f(centr,p)
      neval = 21
      ier = 1
!
!  compute the integral using the 10- and 21-point formula.
!
      do 70 l = 1,3
      go to (5,25,45),l
    5 res10 = C_0_R
      res21 = w21b(6)*fcentr
      resabs = w21b(6)*abs(fcentr)
      do 10 k=1,5
        absc = hlgth*x1(k)
        fval1 = f(centr+absc,p)
        fval2 = f(centr-absc,p)
        fval = fval1+fval2
        res10 = res10+w10(k)*fval
        res21 = res21+w21a(k)*fval
        resabs = resabs+w21a(k)*(abs(fval1)+abs(fval2))
        savfun(k) = fval
        fv1(k) = fval1
        fv2(k) = fval2
   10 continue
      ipx = 5
      do 15 k=1,5
        ipx = ipx+1
        absc = hlgth*x2(k)
        fval1 = f(centr+absc,p)
        fval2 = f(centr-absc,p)
        fval = fval1+fval2
        res21 = res21+w21b(k)*fval
        resabs = resabs+w21b(k)*(abs(fval1)+abs(fval2))
        savfun(ipx) = fval
        fv3(k) = fval1
        fv4(k) = fval2
   15 continue
!
!  test for convergence.
!
      result = res21*hlgth
      resabs = resabs*dhlgth
      reskh = C_05_R*res21
      resasc = w21b(6)*abs(fcentr-reskh)
      do 20 k = 1,5
        resasc = resasc+w21a(k)*(abs(fv1(k)-reskh)+abs(fv2(k)-reskh))&
                        +w21b(k)*(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
   20 continue
      abserr = abs((res21-res10)*hlgth)
      resasc = resasc*dhlgth
      go to 65
!
!  compute the integral using the 43-point formula.
!
   25 res43 = w43b(12)*fcentr
      neval = 43
      do 30 k=1,10
        res43 = res43+savfun(k)*w43a(k)
   30 continue
      do 40 k=1,11
        ipx = ipx+1
        absc = hlgth*x3(k)
        fval = f(absc+centr,p)+f(centr-absc,p)
        res43 = res43+fval*w43b(k)
        savfun(ipx) = fval
   40 continue
!
!  test for convergence.
!
      result = res43*hlgth
      abserr = abs((res43-res21)*hlgth)
      go to 65
!
!  compute the integral using the 87-point formula.
!
   45 res87 = w87b(23)*fcentr
      neval = 87
      do 50 k=1,21
        res87 = res87+savfun(k)*w87a(k)
   50 continue
      do 60 k=1,22
        absc = hlgth*x4(k)
        res87 = res87+w87b(k)*(f(absc+centr,p)+f(centr-absc,p))
   60 continue
      result = res87*hlgth
      abserr = abs((res87-res43)*hlgth)
   65 if(resasc.ne.C_0_R.and.abserr.ne. C_0_R)&
        abserr = resasc*min(C_1_R,&
        (C_200_L*abserr/resasc)**C_1P5_L)
      if (resabs.gt.uflow/(C_50_L*epmach)) abserr = MAX&
        ((epmach*C_50_L)*resabs,abserr)
      if (abserr.le.MAX(epsabs,epsrel*abs(result))) ier = 0
! ***jump out of do-loop
      if (ier.eq.0) go to 999
   70 continue
   80 call xerror(26habnormal return from  qng ,26,ier,0)
  999 return
      end
