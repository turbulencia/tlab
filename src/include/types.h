
#ifndef TYPES_H_INCLUDE
#define TYPES_H_INCLUDE

#define SIZEOFREAL 8
#define SIZEOFINT  4
#define SIZEOFLONGINT  8

#define TREAL      REAL(SIZEOFREAL)
#define TINTEGER   INTEGER(SIZEOFINT)
#define TCOMPLEX   COMPLEX(SIZEOFREAL)
#define TLONGINTEGER   INTEGER(SIZEOFLONGINT)

#define G_FORMAT_R  E13.5E3

#define USE_ACCESS_STREAM

#ifdef PARALLEL

#define GEN_ONEOUTPUT

#ifdef CRAY
#undef USE_ACCESS_STREAM
#define SINGLE_PREC
#define USE_RECLEN
#endif

#ifdef NEC
#undef USE_ACCESS_STREAM
#endif

#ifdef HITACHI
#undef USE_ACCESS_STREAM
#define USE_BUFLESS
#endif

#ifdef SUN
#undef USE_ACCESS_STREAM
#define USE_BINARY
#endif

#endif


#ifdef SINGLE_PREC

#define M_MACH     r1mach
#define M_REAL(x)  FLOAT(x)
#define C_1EM13_R  1.0e-13
#define C_1EM10_R  1.0e-10
#define C_1EM8_R   1.0e-8
#define C_12EM7_R  1.2e-7
#define C_1EM6_R   1.0e-6
#define C_1EM5_R   1.0e-5
#define C_1EM4_R   1.0e-4
#define C_1EM3_R   1.0e-3
#define C_1EM2_R   1.0e-2
#define C_1EM1_R   1.0e-1
#define C_PI_R     3.14159265358979323846
#define C_005_R    5.0e-2
#define C_01_R     0.1
#define C_0125_R   0.125
#define C_025_R    0.25
#define C_05_R     0.5
#define C_075_R    0.75
#define C_0_R      0.0
#define C_1_R      1.0
#define C_1_5_R    1.5
#define C_2_R      2.0
#define C_3_R      3.0
#define C_4_R      4.0
#define C_5_R      5.0
#define C_6_R      6.0
#define C_7_R      7.0
#define C_8_R      8.0
#define C_9_R      9.0

#define C_10_R     10.0
#define C_11_R     11.0
#define C_12_R     12.0
#define C_13_R     13.0
#define C_14_R     14.0
#define C_15_R     15.0
#define C_16_R     16.0
#define C_17_R     17.0
#define C_18_R     18.0
#define C_20_R     20.0
#define C_22_R     22.0
#define C_25_R     25.0
#define C_27_R     27.0
#define C_32_R     32.0
#define C_36_R     36.0
#define C_44_R     44.0
#define C_50_R     50.0
#define C_51_R     51.0
#define C_56_R     56.0
#define C_62_R     62.0

#define C_100_R    100.0
#define C_101_R    101.0
#define C_111_R    111.0
#define C_128_R    128.0
#define C_145_R    145.0
#define C_150_R    150.0
#define C_153_R    153.0
#define C_256_R    256.0

#define C_1E1_R    1.0E+1
#define C_1E2_R    1.0E+2
#define C_1E3_R    1.0E+3
#define C_1E4_R    1.0E+4
#define C_1E5_R    1.0E+5
#define C_1E6_R    1.0E+6

#define C_1432997174477_R  1432997174477.0
#define C_9575080441755_R  9575080441755.0
#define C_5161836677717_R  5161836677717.0
#define C_13612068292357_R 13612068292357.0
#define C_1720146321549_R  1720146321549.0
#define C_2090206949498_R  2090206949498.0
#define C_3134564353537_R  3134564353537.0
#define C_4481467310338_R  4481467310338.0
#define C_2277821191437_R  2277821191437.0
#define C_14882151754819_R 14882151754819.0
			   
#define C_2526269341429_R  2526269341429.0
#define C_6820363962896_R  6820363962896.0
#define C_2006345519317_R  2006345519317.0
#define C_3224310063776_R  3224310063776.0
#define C_2802321613138_R  2802321613138.0
#define C_2924317926251_R  2924317926251.0
               	           
#define C_567301805773_R   567301805773.0
#define C_1357537059087_R  1357537059087.0
#define C_2404267990393_R  2404267990393.0
#define C_2016746695238_R  2016746695238.0
#define C_3550918686646_R  3550918686646.0
#define C_2091501179385_R  2091501179385.0
#define C_1275806237668_R  1275806237668.0
#define C_842570457699_R   842570457699.0

#define C_SMALL_R          1.0e-20
#define C_BIG_R            1.0e20

#define C_OE0_R             0.0745674
#define C_OP0_R             0.222841
#define C_OP1_R             0.0966586
#define C_OP2_R             0.080679
#define C_OP3_R             0.0425497
#define C_OP4_R             0.0219527
#define C_OP5_R             0.00537947

#define C_BE0_R             0.00982802
#define C_BP0_R             0.197301
#define C_BP1_R             0.207126
#define C_BP2_R             0.048025
#define C_BP3_R             0.067662
#define C_BP4_R             0.0200048
#define C_BP5_R             0.0015502
#define C_BR0_R             0.338223

#else

#define M_MACH     d1mach
#define M_REAL(x)  DBLE(x)
#define C_1EM13_R  1.0d-13
#define C_1EM10_R  1.0d-10
#define C_1EM8_R   1.0d-8
#define C_12EM7_R  1.2d-7
#define C_1EM6_R   1.0d-6
#define C_1EM5_R   1.0d-5
#define C_1EM4_R   1.0d-4
#define C_1EM3_R   1.0d-3
#define C_1EM2_R   1.0d-2
#define C_1EM1_R   1.0d-1
#define C_PI_R     3.14159265358979323846d0
#define C_005_R    5.0d-2
#define C_01_R     0.1d0
#define C_0125_R   0.125d0
#define C_025_R    0.25d0
#define C_05_R     0.5d0
#define C_075_R    0.75d0
#define C_0_R      0.0d0
#define C_1_R      1.0d0
#define C_1_5_R    1.5d0
#define C_2_R      2.0d0
#define C_3_R      3.0d0
#define C_4_R      4.0d0
#define C_5_R      5.0d0
#define C_6_R      6.0d0
#define C_7_R      7.0d0
#define C_8_R      8.0d0
#define C_9_R      9.0d0
		   
#define C_10_R     10.0d0
#define C_11_R     11.0d0
#define C_12_R     12.0d0
#define C_13_R     13.0d0
#define C_14_R     14.0d0
#define C_15_R     15.0d0
#define C_16_R     16.0d0
#define C_17_R     17.0d0
#define C_18_R     18.0d0
#define C_20_R     20.0d0
#define C_22_R     22.0d0
#define C_25_R     25.0d0
#define C_27_R     27.0d0
#define C_32_R     32.0d0
#define C_36_R     36.0d0
#define C_44_R     44.0d0
#define C_50_R     50.0d0
#define C_51_R     51.0d0
#define C_56_R     56.0d0
#define C_62_R     62.0d0

#define C_100_R    100.0d0
#define C_101_R    101.0d0
#define C_111_R    111.0d0
#define C_128_R    128.0d0
#define C_145_R    145.0d0
#define C_150_R    150.0d0
#define C_153_R    153.0d0
#define C_256_R    256.0d0

#define C_1E1_R    1.0D+1
#define C_1E2_R    1.0D+2
#define C_1E3_R    1.0D+3
#define C_1E4_R    1.0D+4
#define C_1E5_R    1.0D+5
#define C_1E6_R    1.0D+6

#define C_1432997174477_R  1432997174477.0d0
#define C_9575080441755_R  9575080441755.0d0
#define C_5161836677717_R  5161836677717.0d0
#define C_13612068292357_R 13612068292357.0d0
#define C_1720146321549_R  1720146321549.0d0
#define C_2090206949498_R  2090206949498.0d0
#define C_3134564353537_R  3134564353537.0d0
#define C_4481467310338_R  4481467310338.0d0
#define C_2277821191437_R  2277821191437.0d0
#define C_14882151754819_R 14882151754819.0d0
			   
#define C_2526269341429_R  2526269341429.0d0
#define C_6820363962896_R  6820363962896.0d0
#define C_2006345519317_R  2006345519317.0d0
#define C_3224310063776_R  3224310063776.0d0
#define C_2802321613138_R  2802321613138.0d0
#define C_2924317926251_R  2924317926251.0d0
               	           
#define C_567301805773_R   567301805773.0d0
#define C_1357537059087_R  1357537059087.0d0
#define C_2404267990393_R  2404267990393.0d0
#define C_2016746695238_R  2016746695238.0d0
#define C_3550918686646_R  3550918686646.0d0
#define C_2091501179385_R  2091501179385.0d0
#define C_1275806237668_R  1275806237668.0d0
#define C_842570457699_R   842570457699.0d0

#define C_SMALL_R          1.0d-20
#define C_BIG_R            1.0d20

#define C_OE0_R             0.0745674d0
#define C_OP0_R             0.222841d0
#define C_OP1_R             0.0966586d0
#define C_OP2_R             0.080679d0
#define C_OP3_R             0.0425497d0
#define C_OP4_R             0.0219527d0
#define C_OP5_R             0.00537947d0

#define C_BE0_R             0.00982802d0
#define C_BP0_R             0.197301d0
#define C_BP1_R             0.207126d0
#define C_BP2_R             0.048025d0
#define C_BP3_R             0.067662d0
#define C_BP4_R             0.0200048d0
#define C_BP5_R             0.0015502d0
#define C_BR0_R             0.338223d0

#endif

#endif			   
