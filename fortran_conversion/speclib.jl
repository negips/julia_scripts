#using FortranFiles
using OffsetArrays
using Parameters
using Printf

# =============================================================================
#
#     LIBRARY ROUTINES FOR SPECTRAL METHODS
#
#     March 1989
#
#     For questions, comments or suggestions, please contact:
#
#     Einar Malvin Ronquist
#     Room 3-243
#     Department of Mechanical Engineering
#     Massachusetts Institute of Technology
#     77 Massachusetts Avenue
#     Cambridge, MA 0299
#     U.S.A.
#
#------------------------------------------------------------------------------
#
#     ABBRIVIATIONS:
#
#     M   - Set of mesh points
#     Z   - Set of collocation/quadrature points
#     W   - Set of quadrature weights
#     H   - Lagrangian interpolant
#     D   - Derivative operator
#     I   - Interpolation operator
#     GL  - Gauss Legendre
#     GLL - Gauss-Lobatto Legendre
#     GJ  - Gauss Jacobi
#     GLJ - Gauss-Lobatto Jacobi
#
#
#     MAIN ROUTINES:
#
#     Points and weights:
#
#     ZWGL      Compute Gauss Legendre points and weights
#     ZWGLL     Compute Gauss-Lobatto Legendre points and weights
#     ZWGJ      Compute Gauss Jacobi points and weights (general)
#     ZWGLJ     Compute Gauss-Lobatto Jacobi points and weights (general)
#
#     Lagrangian interpolants:
#
#     HGL       Compute Gauss Legendre Lagrangian interpolant
#     HGLL      Compute Gauss-Lobatto Legendre Lagrangian interpolant
#     HGJ       Compute Gauss Jacobi Lagrangian interpolant (general)
#     HGLJ      Compute Gauss-Lobatto Jacobi Lagrangian interpolant (general)
#
#     Derivative operators:
#
#     DGLL      Compute Gauss-Lobatto Legendre derivative matrix
#     DGLLGL    Compute derivative matrix for a staggered mesh (GLL->GL)
#     DGJ       Compute Gauss Jacobi derivative matrix (general)
#     DGLJ      Compute Gauss-Lobatto Jacobi derivative matrix (general)
#     DGLJGJ    Compute derivative matrix for a staggered mesh (GLJ->GJ) (general)
#
#     Interpolation operators:
#
#     IGLM      Compute interpolation operator GL  -> M
#     IGLLM     Compute interpolation operator GLL -> M
#     IGJM      Compute interpolation operator GJ  -> M  (general)
#     IGLJM     Compute interpolation operator GLJ -> M  (general)
#
#     Other:
#
#     PNLEG     Compute Legendre polynomial of degree N
#     PNDLEG    Compute derivative of Legendre polynomial of degree N
#
#     Comments:
#
#     Note that many of the above routines exist in both single and
#     double precision. If the name of the single precision routine is
#     SUB, the double precision version is called SUBD. In most cases
#     all the "low-level" arithmetic is done in double precision, even
#     for the single precsion versions.
#
#     Useful references:
#
# [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
#     Providence, Rhode Island, 1939.
# [2] Abramowitz & Stegun: Handbook of Mathematical Functions,
#     Dover, New York, 1972.
# [3] Canuto, Hussaini, Quarteroni & Zang: Spectral Methods in Fluid
#     Dynamics, Springer-Verlag, 1988.
#
#
# =============================================================================
#
#--------------------------------------------------------------------
      function ZWGL(Z, W, NP)
#--------------------------------------------------------------------
#
#     Generate NP Gauss Legendre points (Z) and weights (W)
#     associated with Jacobi polynomial P(N)(alpha=0,beta=0).
#     The polynomial degree N=NP-1.
#     Z and W are in single precision, but all the arithmetic
#     operations are done in double precision.
#
#--------------------------------------------------------------------
#      REAL Z[1],W[1]
      ALPHA = 0.
      BETA  = 0.
      ZWGJ(Z, W, NP, ALPHA, BETA)
      return nothing
      end
#
      function ZWGLL(Z, W, NP)
#--------------------------------------------------------------------
#
#     Generate NP Gauss-Lobatto Legendre points (Z) and weights (W)
#     associated with Jacobi polynomial P(N)(alpha=0,beta=0).
#     The polynomial degree N=NP-1.
#     Z and W are in single precision, but all the arithmetic
#     operations are done in double precision.
#
#--------------------------------------------------------------------
#      REAL Z[1],W[1]
      ALPHA = 0.
      BETA  = 0.
      ZWGLJ(Z, W, NP, ALPHA, BETA)
      return nothing
      end
#
      function ZWGJ(Z, W, NP, ALPHA, BETA)
#--------------------------------------------------------------------
#
#     Generate NP GAUSS JACOBI points (Z) and weights (W)
#     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
#     The polynomial degree N=NP-1.
#     Single precision version.
#
#--------------------------------------------------------------------
                 NMAX = 84
                 lzd = NMAX
       
                 ZD = zeros(typeof(Z[1]),NMAX)
#      REAL*8  ZD[lzd],WD[lzd],APHAD,BETAD
#      REAL Z[1],W[1],ALPHA,BETA
#
      NPMAX = lzd
      if NP > NPMAX
         println(stdout, "Too large polynomial degree in ZWGJ")
         println(stdout, "Maximum polynomial degree is", NMAX)
         println(stdout, "Here NP=", NP)
         exitt
      end
      ALPHAD = ALPHA
      BETAD  = BETA
      ZWGJD(ZD, WD, NP, ALPHAD, BETAD)
      for I = 1:NP
         Z[I] = ZD[I]
         W[I] = WD[I]
      end
      return nothing
      end
#
      function ZWGJD(Z, W, NP, ALPHA, BETA)
#--------------------------------------------------------------------
#
#     Generate NP GAUSS JACOBI points (Z) and weights (W)
#     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
#     The polynomial degree N=NP-1.
#     Double precision version.
#
#--------------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  Z[1],W[1],ALPHA,BETA
#
      N     = NP-1
      DN    = ((N))
      ONE   = 1.
      TWO   = 2.
      APB   = ALPHA+BETA
#
      if NP <= 0
         println(stdout, "ZWGJD: Minimum number of Gauss points is 1", np)
         exitt
      end
      if (ALPHA <= -ONE) || (BETA <= -ONE)
         println(stdout, "ZWGJD: Alpha and Beta must be greater than -1")
         exitt
      end
#
      if NP == 1
         Z[1] = (BETA-ALPHA)/(APB+TWO)
         W[1] = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)/GAMMAF(APB+TWO)*
                 TWO^(APB+ONE)
         return nothing
      end
#
      JACG(Z, NP, ALPHA, BETA)
#
      NP1   = N+1
      NP2   = N+2
      DNP1  = ((NP1))
      DNP2  = ((NP2))
      FAC1  = DNP1+ALPHA+BETA+ONE
      FAC2  = FAC1+DNP1
      FAC3  = FAC2+ONE
      FNORM = PNORMJ(NP1, ALPHA, BETA)
      RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)
      for I = 1:NP
         JACOBF(P, PD, PM1, PDM1, PM2, PDM2, NP2, ALPHA, BETA, Z[I])
         W[I] = -RCOEF/(P*PDM1)
      end
      return nothing
      end
#
      function ZWGLJ(Z, W, NP, ALPHA, BETA)
#--------------------------------------------------------------------
#
#     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
#     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
#     The polynomial degree N=NP-1.
#     Single precision version.
#
#--------------------------------------------------------------------
                 NMAX = 84
                 lzd = NMAX
                 ZD = zeros(typeof(Z[1]),NMAX)
                 WD = zeros(typeof(W[1]),NMAX)
              
#      REAL*8  ZD[lzd],WD[lzd],ALPHAD,BETAD
#      REAL Z[1],W[1],ALPHA,BETA
#
      NPMAX = lzd
      if NP > NPMAX
         println(stdout, "Too large polynomial degree in ZWGLJ")
         println(stdout, "Maximum polynomial degree is", NMAX)
         println(stdout, "Here NP=", NP)
         exitt
      end
      ALPHAD = ALPHA
      BETAD  = BETA
      ZWGLJD(ZD, WD, NP, ALPHAD, BETAD)
      for I = 1:NP
         Z[I] = ZD[I]
         W[I] = WD[I]
      end
      return nothing
      end
#
      function ZWGLJD(Z, W, NP, ALPHA, BETA)
#--------------------------------------------------------------------
#
#     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
#     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
#     The polynomial degree N=NP-1.
#     Double precision version.
#
#--------------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  Z[NP],W[NP],ALPHA,BETA
#
      N     = NP-1
      NM1   = N-1
      ONE   = 1.
      TWO   = 2.
#
      if NP <= 1
       println(stdout, "ZWGLJD: Minimum number of Gauss-Lobatto points is 2")
       println(stdout, "ZWGLJD: alpha,beta:", alpha, beta, np)
       exitt
      end
      if (ALPHA <= -ONE) || (BETA <= -ONE)
         println(stdout, "ZWGLJD: Alpha and Beta must be greater than -1")
         exitt
      end
#
      if NM1 > 0
         ALPG  = ALPHA+ONE
         BETG  = BETA+ONE
         ZWGJD(Z[2], W[2], NM1, ALPG, BETG)
      end
      Z[1]  = -ONE
      Z[NP] =  ONE
      for I = 2:NP-1
         W[I] = W[I]/(ONE-Z[I]^2)
      end
      JACOBF(P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z[1])
      W[1]  = ENDW1(N, ALPHA, BETA)/(TWO*PD)
      JACOBF(P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z[NP])
      W[NP] = ENDW2(N, ALPHA, BETA)/(TWO*PD)
#
      return nothing
      end
#
      function ENDW1(N, ALPHA, BETA)
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  ALPHA,BETA
      ZERO  = 0.
      ONE   = 1.
      TWO   = 2.
      THREE = 3.
      FOUR  = 4.
      APB   = ALPHA+BETA
      if N == 0
         ENDW1 = ZERO
         return nothing
      end
      F1   = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO^(APB+TWO)/TWO
      if N == 1
         ENDW1 = F1
         return nothing
      end
      FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO^(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO^(APB+THREE)
      F2    = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2)*
               (APB+THREE)/FOUR
      if N == 2
         ENDW1 = F2
         return nothing
      end
      for I = 3:N
         DI   = ((I-1))
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
      end
      ENDW1  = F3
      return nothing
      end
#
      function ENDW2(N, ALPHA, BETA)
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  ALPHA,BETA
      ZERO  = 0.
      ONE   = 1.
      TWO   = 2.
      THREE = 3.
      FOUR  = 4.
      APB   = ALPHA+BETA
      if N == 0
         ENDW2 = ZERO
         return nothing
      end
      F1   = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO^(APB+TWO)/TWO
      if N == 1
         ENDW2 = F1
         return nothing
      end
      FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO^(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO^(APB+THREE)
      F2    = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2)*
               (APB+THREE)/FOUR
      if N == 2
         ENDW2 = F2
         return nothing
      end
      for I = 3:N
         DI   = ((I-1))
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   =  -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
      end
      ENDW2  = F3
      return nothing
      end
#
      function GAMMAF(X)
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  X
      ZERO = 0.0
      HALF = 0.5
      ONE  = 1.0
      TWO  = 2.0
      FOUR = 4.0
      PI   = FOUR*atan(ONE)
      GAMMAF = ONE
      if (X == -HALF) GAMMAF = -TWO*sqrt(PI) end
      if (X == HALF) GAMMAF =  sqrt(PI) end
      if (X == ONE) GAMMAF =  ONE end
      if (X == TWO) GAMMAF =  ONE end
      if (X == 1.5) GAMMAF =  sqrt(PI)/2. end
      if (X == 2.5) GAMMAF =  1.5*sqrt(PI)/2. end
      if (X == 3.5) GAMMAF =  0.5*(2.5*(1.5*sqrt(PI))) end
      if (X == 3.) GAMMAF =  2. end
      if (X == 4.) GAMMAF = 6. end
      if (X == 5.) GAMMAF = 24. end
      if (X == 6.) GAMMAF = 120. end
      return nothing
      end
#
      function PNORMJ(N, ALPHA, BETA)
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  ALPHA,BETA
      ONE   = 1.
      TWO   = 2.
      DN    = ((N))
      CONST = ALPHA+BETA+ONE
      if N <= 1
         PROD   = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA)
         PROD   = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA))
         PNORMJ = PROD * TWO^CONST/(TWO*DN+CONST)
         return nothing
      end
      PROD  = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)
      PROD  = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE))
      PROD  = PROD*(ONE+ALPHA)*(TWO+ALPHA)
      PROD  = PROD*(ONE+BETA)*(TWO+BETA)
      for I = 3:N
         DINDX = ((I))
         FRAC  = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA))
         PROD  = PROD*FRAC
      end
      PNORMJ = PROD * TWO^CONST/(TWO*DN+CONST)
      return nothing
      end
#
      function JACG(XJAC, NP, ALPHA, BETA)
#--------------------------------------------------------------------
#
#     Compute NP Gauss points XJAC, which are the zeros of the
#     Jacobi polynomial J(NP) with parameters ALPHA and BETA.
#     ALPHA and BETA determines the specific type of Gauss points.
#     Examples:
#     ALPHA = BETA =  0.0  ->  Legendre points
#     ALPHA = BETA = -0.5  ->  Chebyshev points
#
#--------------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  XJAC[1]
           KSTOP  = 10
           EPS = 1.0f-12
      N   = NP-1
      ONE = 1.
      DTH = 4.0*atan(ONE)/(2.0*((N))+2.)
      for J = 1:NP
         if J == 1
            X = cos((2.0*(((J))-1.)+1.)*DTH)
         else
            X1 = cos((2.0*(((J))-1.)+1.)*DTH)
            X2 = XLAST
            X  = (X1+X2)/2.
         end
         for K = 1:KSTOP
            JACOBF(P, PD, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, X)
            RECSUM = 0.
            JM = J-1
            for I = 1:JM
               RECSUM = RECSUM+1.0/(X-XJAC[NP-I+1])
            end
            DELX = -P/(PD-RECSUM*P)
            X    = X+DELX
            if (abs(DELX) < EPS) @goto L31 end
         end
         @label L31
         XJAC[NP-J+1] = X
         XLAST        = X
      end
      for I = 1:NP
         XMIN = 2.
         for J = I:NP
            if XJAC[J] < XMIN
               XMIN = XJAC[J]
               JMIN = J
            end
         end
         if JMIN != I
            SWAP = XJAC[I]
            XJAC[I] = XJAC[JMIN]
            XJAC[JMIN] = SWAP
         end
      end
      return nothing
      end
#
      function JACOBF(POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2,
                         N, ALP, BET, X)
#--------------------------------------------------------------------
#
#     Computes the Jacobi polynomial (POLY) and its derivative (PDER)
#     of degree N at X.
#
#--------------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
      APB  = ALP+BET
      POLY = 1.
      PDER = 0.
      if (N == 0) return end
      POLYL = POLY
      PDERL = PDER
      POLY  = (ALP-BET+(APB+2.)*X)/2.
      PDER  = (APB+2.)/2.
      if (N == 1) return end
      for K = 2:N
         DK = ((K))
         A1 = 2.0*DK*(DK+APB)*(2.0*DK+APB-2.)
         A2 = (2.0*DK+APB-1.)*(ALP^2-BET^2)
         B3 = (2.0*DK+APB-2.)
         A3 = B3*(B3+1.)*(B3+2.)
         A4 = 2.0*(DK+ALP-1.)*(DK+BET-1.)*(2.0*DK+APB)
         POLYN  = ((A2+A3*X)*POLY-A4*POLYL)/A1
         PDERN  = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1
         PSAVE  = POLYL
         PDSAVE = PDERL
         POLYL  = POLY
         POLY   = POLYN
         PDERL  = PDER
         PDER   = PDERN
      end
      POLYM1 = POLYL
      PDERM1 = PDERL
      POLYM2 = PSAVE
      PDERM2 = PDSAVE
      return nothing
      end
#
      function HGJ(II, Z, ZGJ, NP, ALPHA, BETA)
#---------------------------------------------------------------------
#
#     Compute the value of the Lagrangian interpolant HGJ through
#     the NP Gauss Jacobi points ZGJ at the point Z.
#     Single precision version.
#
#---------------------------------------------------------------------
                 NMAX = 84
                 lzd = NMAX
#      REAL*8  ZD,ZGJD[lzd],HGJD,ALPHAD,BETAD
#      REAL Z,ZGJ[1],ALPHA,BETA
      NPMAX = lzd
      if NP > NPMAX
         println(stdout, "Too large polynomial degree in HGJ")
         println(stdout, "Maximum polynomial degree is", NMAX)
         println(stdout, "Here NP=", NP)
         exitt
      end
      ZD = Z
      for I = 1:NP
         ZGJD[I] = ZGJ[I]
      end
      ALPHAD = ALPHA
      BETAD  = BETA
      HGJ    = HGJD(II, ZD, ZGJD, NP, ALPHAD, BETAD)
      return nothing
      end
#
      function HGJD(II, Z, ZGJ, NP, ALPHA, BETA)
#---------------------------------------------------------------------
#
#     Compute the value of the Lagrangian interpolant HGJD through
#     the NZ Gauss-Lobatto Jacobi points ZGJ at the point Z.
#     Double precision version.
#
#---------------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  Z,ZGJ[1],ALPHA,BETA
      EPS = 1.0f-5
      ONE = 1.
      ZI  = ZGJ[II]
      DZ  = Z-ZI
      if abs(DZ) < EPS
         HGJD = ONE
         return nothing
      end
      JACOBF(PZI, PDZI, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, ZI)
      JACOBF(PZ, PDZ, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, Z)
      HGJD  = PZ/(PDZI*(Z-ZI))
      return nothing
      end
#
      function HGLJ(II, Z, ZGLJ, NP, ALPHA, BETA)
#---------------------------------------------------------------------
#
#     Compute the value of the Lagrangian interpolant HGLJ through
#     the NZ Gauss-Lobatto Jacobi points ZGLJ at the point Z.
#     Single precision version.
#
#---------------------------------------------------------------------
                 NMAX = 84
                 lzd = NMAX
#      REAL*8  ZD,ZGLJD[lzd],HGLJD,ALPHAD,BETAD
#      REAL Z,ZGLJ[1],ALPHA,BETA
      NPMAX = lzd
      if NP > NPMAX
         println(stdout, "Too large polynomial degree in HGLJ")
         println(stdout, "Maximum polynomial degree is", NMAX)
         println(stdout, "Here NP=", NP)
         exitt
      end
      ZD = Z
      for I = 1:NP
         ZGLJD[I] = ZGLJ[I]
      end
      ALPHAD = ALPHA
      BETAD  = BETA
      HGLJ   = HGLJD(II, ZD, ZGLJD, NP, ALPHAD, BETAD)
      return nothing
      end
#
      function HGLJD(I, Z, ZGLJ, NP, ALPHA, BETA)
#---------------------------------------------------------------------
#
#     Compute the value of the Lagrangian interpolant HGLJD through
#     the NZ Gauss-Lobatto Jacobi points ZJACL at the point Z.
#     Double precision version.
#
#---------------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  Z,ZGLJ[1],ALPHA,BETA
      EPS = 1.0f-5
      ONE = 1.
      ZI  = ZGLJ[I]
      DZ  = Z-ZI
      if abs(DZ) < EPS
         HGLJD = ONE
         return nothing
      end
      N      = NP-1
      DN     = ((N))
      EIGVAL = -DN*(DN+ALPHA+BETA+ONE)
      JACOBF(PI, PDI, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, ZI)
      CONST  = EIGVAL*PI+ALPHA*(ONE+ZI)*PDI-BETA*(ONE-ZI)*PDI
      JACOBF(P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z)
      HGLJD  = (ONE-Z^2)*PD/(CONST*(Z-ZI))
      return nothing
      end
#
      function DGJ(D, DT, Z, NZ, lzd, ALPHA, BETA)
#-----------------------------------------------------------------
#
#     Compute the derivative matrix D and its transpose DT
#     associated with the Nth order Lagrangian interpolants
#     through the NZ Gauss Jacobi points Z.
#     Note: D and DT are square matrices.
#     Single precision version.
#
#-----------------------------------------------------------------
                 NMAX = 84
                 lzdD = NMAX
#      REAL*8  DD[lzdD,lzdD],DTD[lzdD,lzdD],ZD[lzdD],ALPHAD,BETAD
#      REAL D[lzd,lzd],DT[lzd,lzd],Z[1],ALPHA,BETA
#
      if NZ <= 0
         println(stdout, "DGJ: Minimum number of Gauss points is 1")
         exitt
      end
      if NZ > NMAX
         println(stdout, "Too large polynomial degree in DGJ")
         println(stdout, "Maximum polynomial degree is", NMAX)
         println(stdout, "Here Nz=", Nz)
         exitt
      end
      if (ALPHA <= -1.) || (BETA <= -1.)
         println(stdout, "DGJ: Alpha and Beta must be greater than -1")
         exitt
      end
      ALPHAD = ALPHA
      BETAD  = BETA
      for I = 1:NZ
         ZD[I] = Z[I]
      end
      DGJD(DD, DTD, ZD, NZ, lzdD, ALPHAD, BETAD)
      for I = 1:NZ
      for J = 1:NZ
         D[I,J]  = DD[I,J]
         DT[I,J] = DTD[I,J]
      end
      end
      return nothing
      end
#
      function DGJD(D, DT, Z, NZ, lzd, ALPHA, BETA)
#-----------------------------------------------------------------
#
#     Compute the derivative matrix D and its transpose DT
#     associated with the Nth order Lagrangian interpolants
#     through the NZ Gauss Jacobi points Z.
#     Note: D and DT are square matrices.
#     Double precision version.
#
#-----------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  D[lzd,lzd],DT[lzd,lzd],Z[1],ALPHA,BETA
      N    = NZ-1
      DN   = ((N))
      ONE  = 1.
      TWO  = 2.
#
      if NZ <= 1
       println(stdout, "DGJD: Minimum number of Gauss-Lobatto points is 2")
       exitt
      end
      if (ALPHA <= -ONE) || (BETA <= -ONE)
         println(stdout, "DGJD: Alpha and Beta must be greater than -1")
         exitt
      end
#
      for I = 1:NZ
      for J = 1:NZ
         JACOBF(PI, PDI, PM1, PDM1, PM2, PDM2, NZ, ALPHA, BETA, Z[I])
         JACOBF(PJ, PDJ, PM1, PDM1, PM2, PDM2, NZ, ALPHA, BETA, Z[J])
         if (I != J) D[I,J] = PDI/(PDJ*(Z[I]-Z[J])) end
         if (I == J) D[I,J] = ((ALPHA+BETA+TWO)*Z[I]+ALPHA-BETA)/(
                              TWO*(ONE-Z[I]^2))
         end
         DT[J,I] = D[I,J]
      end
      end
      return nothing
      end
#
      function DGLJ(D, DT, Z, NZ, lzd, ALPHA, BETA)
#-----------------------------------------------------------------
#
#     Compute the derivative matrix D and its transpose DT
#     associated with the Nth order Lagrangian interpolants
#     through the NZ Gauss-Lobatto Jacobi points Z.
#     Note: D and DT are square matrices.
#     Single precision version.
#
#-----------------------------------------------------------------
                 NMAX = 84
                 lzdD = NMAX
#      REAL*8  DD[lzdD,lzdD],DTD[lzdD,lzdD],ZD[lzdD],ALPHAD,BETAD
#      REAL D[lzd,lzd],DT[lzd,lzd],Z[1],ALPHA,BETA
#
      if NZ <= 1
       println(stdout, "DGLJ: Minimum number of Gauss-Lobatto points is 2")
       exitt
      end
      if NZ > NMAX
         println(stdout, "Too large polynomial degree in DGLJ")
         println(stdout, "Maximum polynomial degree is", NMAX)
         println(stdout, "Here NZ=", NZ)
         exitt
      end
      if (ALPHA <= -1.) || (BETA <= -1.)
         println(stdout, "DGLJ: Alpha and Beta must be greater than -1")
         exitt
      end
      ALPHAD = ALPHA
      BETAD  = BETA
      for I = 1:NZ
         ZD[I] = Z[I]
      end
      DGLJD(DD, DTD, ZD, NZ, lzdD, ALPHAD, BETAD)
      for I = 1:NZ
      for J = 1:NZ
         D[I,J]  = DD[I,J]
         DT[I,J] = DTD[I,J]
      end
      end
      return nothing
      end
#
      function DGLJD(D, DT, Z, NZ, lzd, ALPHA, BETA)
#-----------------------------------------------------------------
#
#     Compute the derivative matrix D and its transpose DT
#     associated with the Nth order Lagrangian interpolants
#     through the NZ Gauss-Lobatto Jacobi points Z.
#     Note: D and DT are square matrices.
#     Double precision version.
#
#-----------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  D[lzd,lzd],DT[lzd,lzd],Z[1],ALPHA,BETA
      N    = NZ-1
      DN   = ((N))
      ONE  = 1.
      TWO  = 2.
      EIGVAL = -DN*(DN+ALPHA+BETA+ONE)
#
      if NZ <= 1
       println(stdout, "DGLJD: Minimum number of Gauss-Lobatto points is 2")
       exitt
      end
      if (ALPHA <= -ONE) || (BETA <= -ONE)
         println(stdout, "DGLJD: Alpha and Beta must be greater than -1")
         exitt
      end
#
      for I = 1:NZ
      for J = 1:NZ
         JACOBF(PI, PDI, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z[I])
         JACOBF(PJ, PDJ, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z[J])
         CI = EIGVAL*PI-(BETA*(ONE-Z[I])-ALPHA*(ONE+Z[I]))*PDI
         CJ = EIGVAL*PJ-(BETA*(ONE-Z[J])-ALPHA*(ONE+Z[J]))*PDJ
         if (I != J) D[I,J] = CI/(CJ*(Z[I]-Z[J])) end
         if ((I == J) && (I != 1) && (I != NZ))
         D[I,J] = (ALPHA*(ONE+Z[I])-BETA*(ONE-Z[I]))/(
                  TWO*(ONE-Z[I]^2))
         end
         if ((I == J) && (I == 1))
         D[I,J] =  (EIGVAL+ALPHA)/(TWO*(BETA+TWO))
         end
         if ((I == J) && (I == NZ))
         D[I,J] = -(EIGVAL+BETA)/(TWO*(ALPHA+TWO))
         end
         DT[J,I] = D[I,J]
      end
      end
      return nothing
      end
#
      function DGLL(D, DT, Z, NZ, lzd)
#-----------------------------------------------------------------
#
#     Compute the derivative matrix D and its transpose DT
#     associated with the Nth order Lagrangian interpolants
#     through the NZ Gauss-Lobatto Legendre points Z.
#     Note: D and DT are square matrices.
#
#-----------------------------------------------------------------
                 NMAX = 84
#      REAL D[lzd,lzd],DT[lzd,lzd],Z[1]
      N  = NZ-1
      if NZ > NMAX
         println(stdout, "Subroutine DGLL")
         println(stdout, "Maximum polynomial degree =", NMAX)
         println(stdout, "Polynomial degree         =", NZ)
      end
      if NZ == 1
         D[1,1] = 0.
         return nothing
      end
      FN = (N)
      d0 = FN*(FN+1.)/4.
      for I = 1:NZ
      for J = 1:NZ
         D[I,J] = 0.
         if (I != J) D[I,J] = PNLEG(Z[I], N)/(
                              PNLEG(Z[J], N)*(Z[I]-Z[J]))
         end
         if ((I == J) && (I == 1)) D[I,J] = -d0 end
         if ((I == J) && (I == NZ)) D[I,J] =  d0 end
         DT[J,I] = D[I,J]
      end
      end
      return nothing
      end
#
      function HGLL(I, Z, ZGLL, NZ)
#---------------------------------------------------------------------
#
#     Compute the value of the Lagrangian interpolant L through
#     the NZ Gauss-Lobatto Legendre points ZGLL at the point Z.
#
#---------------------------------------------------------------------
#      REAL ZGLL[1]
      EPS = 1.0f-5
      DZ = Z - ZGLL[I]
      if abs(DZ) < EPS
         HGLL = 1.
         return nothing
      end
      N = NZ - 1
      ALFAN = (N)*((N)+1.)
      HGLL = - (1.0-Z*Z)*PNDLEG(Z, N)/(
               ALFAN*PNLEG(ZGLL[I], N)*(Z-ZGLL[I]))
      return nothing
      end
#
      function HGL(I, Z, ZGL, NZ)
#---------------------------------------------------------------------
#
#     Compute the value of the Lagrangian interpolant HGL through
#     the NZ Gauss Legendre points ZGL at the point Z.
#
#---------------------------------------------------------------------
#      REAL ZGL[1]
      EPS = 1.0f-5
      DZ = Z - ZGL[I]
      if abs(DZ) < EPS
         HGL = 1.
         return nothing
      end
      N = NZ-1
      HGL = PNLEG(Z, NZ)/(PNDLEG(ZGL[I], NZ)*(Z-ZGL[I]))
      return nothing
      end
#
      function PNLEG(Z, N)
#---------------------------------------------------------------------
#
#     Compute the value of the Nth order Legendre polynomial at Z.
#     (Simpler than JACOBF)
#     Based on the recursion formula for the Legendre polynomials.
#
#---------------------------------------------------------------------
#
#     This next statement is to overcome the underflow bug in the i860.  
#     It can be removed at a later date.  11 Aug 1990   pff.
#
      if (abs(Z) < 1.0f-25) Z = 0.0 end
#
      P1   = 1.
      if N == 0
         PNLEG = P1
         return nothing
      end
      P2   = Z
      P3   = P2
      for K = 1:N-1
         FK  = (K)
         P3  = ((2.0*FK+1.)*Z*P2 - FK*P1)/(FK+1.)
         P1  = P2
         P2  = P3
      end
      PNLEG = P3
      if (n == 0) pnleg = 1. end
      return nothing
      end
#
      function PNDLEG(Z, N)
#----------------------------------------------------------------------
#
#     Compute the derivative of the Nth order Legendre polynomial at Z.
#     (Simpler than JACOBF)
#     Based on the recursion formula for the Legendre polynomials.
#
#----------------------------------------------------------------------
      P1   = 1.
      P2   = Z
      P1D  = 0.
      P2D  = 1.
      P3D  = 1.
      for K = 1:N-1
         FK  = (K)
         P3  = ((2.0*FK+1.)*Z*P2 - FK*P1)/(FK+1.)
         P3D = ((2.0*FK+1.)*P2 + (2.0*FK+1.)*Z*P2D - FK*P1D)/(FK+1.)
         P1  = P2
         P2  = P3
         P1D = P2D
         P2D = P3D
      end
      PNDLEG = P3D
      if (N == 0) pndleg = 0. end
      return nothing
      end
#
      function DGLLGL(D, DT, ZM1, ZM2, IM12, NZM1, NZM2, ND1, ND2)
#-----------------------------------------------------------------------
#
#     Compute the (one-dimensional) derivative matrix D and its
#     transpose DT associated with taking the derivative of a variable
#     expanded on a Gauss-Lobatto Legendre mesh (M1), and evaluate its
#     derivative on a Guass Legendre mesh (M2).
#     Need the one-dimensional interpolation operator IM12
#     (see subroutine IGLLGL).
#     Note: D and DT are rectangular matrices.
#
#-----------------------------------------------------------------------
#      REAL D[ND2,ND1], DT[ND1,ND2], ZM1[ND1], ZM2[ND2], IM12[ND2,ND1]
      if NZM1 == 1
        D[1,1] = 0.
        DT[1,1] = 0.
        return nothing
      end
      EPS = 1.0f-6
      NM1 = NZM1-1
      for IP = 1:NZM2
         for JQ = 1:NZM1
            ZP = ZM2[IP]
            ZQ = ZM1[JQ]
            if (abs(ZP) < EPS) && (abs(ZQ) < EPS)
                D[IP,JQ] = 0.
            else
                D[IP,JQ] = (PNLEG(ZP, NM1)/PNLEG(ZQ, NM1)-
                           IM12[IP,JQ])/(ZP-ZQ)
            end
            DT[JQ,IP] = D[IP,JQ]
      end
      end
      return nothing
      end
#
      function DGLJGJ(D, DT, ZGL, ZG, IGLG, NPGL, NPG, ND1, ND2, ALPHA, BETA)
#-----------------------------------------------------------------------
#
#     Compute the (one-dimensional) derivative matrix D and its
#     transpose DT associated with taking the derivative of a variable
#     expanded on a Gauss-Lobatto Jacobi mesh (M1), and evaluate its
#     derivative on a Guass Jacobi mesh (M2).
#     Need the one-dimensional interpolation operator IM12
#     (see subroutine IGLJGJ).
#     Note: D and DT are rectangular matrices.
#     Single precision version.
#
#-----------------------------------------------------------------------
#      REAL D[ND2,ND1], DT[ND1,ND2], ZGL[ND1], ZG[ND2], IGLG[ND2,ND1]
                 NMAX = 84
                 NDD = NMAX
#      REAL*8  DD[NDD,NDD], DTD[NDD,NDD]
#      REAL*8  ZGD[NDD], ZGLD[NDD], IGLGD[NDD,NDD]
#      REAL*8  ALPHAD, BETAD
#
      if NPGL <= 1
       println(stdout, "DGLJGJ: Minimum number of Gauss-Lobatto points is 2")
       exitt
      end
      if NPGL > NMAX
         println(stdout, "Polynomial degree too high in DGLJGJ")
         println(stdout, "Maximum polynomial degree is", NMAX)
         println(stdout, "Here NPGL=", NPGL)
         exitt
      end
      if (ALPHA <= -1.) || (BETA <= -1.)
         println(stdout, "DGLJGJ: Alpha and Beta must be greater than -1")
         exitt
      end
#
      ALPHAD = ALPHA
      BETAD  = BETA
      for I = 1:NPG
         ZGD[I] = ZG[I]
         for J = 1:NPGL
            IGLGD[I,J] = IGLG[I,J]
      end
      end
      for I = 1:NPGL
         ZGLD[I] = ZGL[I]
      end
      DGLJGJD(DD, DTD, ZGLD, ZGD, IGLGD, NPGL, NPG, NDD, NDD, ALPHAD, BETAD)
      for I = 1:NPG
      for J = 1:NPGL
         D[I,J]  = DD[I,J]
         DT[J,I] = DTD[J,I]
      end
      end
      return nothing
      end
#
      function DGLJGJD(D, DT, ZGL, ZG, IGLG, NPGL, NPG, ND1, ND2, ALPHA, BETA)
#-----------------------------------------------------------------------
#
#     Compute the (one-dimensional) derivative matrix D and its
#     transpose DT associated with taking the derivative of a variable
#     expanded on a Gauss-Lobatto Jacobi mesh (M1), and evaluate its
#     derivative on a Guass Jacobi mesh (M2).
#     Need the one-dimensional interpolation operator IM12
#     (see subroutine IGLJGJ).
#     Note: D and DT are rectangular matrices.
#     Double precision version.
#
#-----------------------------------------------------------------------
#      IMPLICIT REAL*8(A-H,O-Z)
#      REAL*8  D[ND2,ND1], DT[ND1,ND2], ZGL[ND1], ZG[ND2]
#      REAL*8  IGLG[ND2,ND1], ALPHA, BETA
#
      if NPGL <= 1
       println(stdout, "DGLJGJD: Minimum number of Gauss-Lobatto points is 2")
       exitt
      end
      if (ALPHA <= -1.) || (BETA <= -1.)
         println(stdout, "DGLJGJD: Alpha and Beta must be greater than -1")
         exitt
      end
#
      EPS    = 1.0f-6
      ONE    = 1.
      TWO    = 2.
      NGL    = NPGL-1
      DN     = ((NGL))
      EIGVAL = -DN*(DN+ALPHA+BETA+ONE)
#
      for I = 1:NPG
      for J = 1:NPGL
         DZ = abs(ZG[I]-ZGL[J])
         if DZ < EPS
            D[I,J] = (ALPHA*(ONE+ZG[I])-BETA*(ONE-ZG[I]))/(
                     TWO*(ONE-ZG[I]^2))
         else
            JACOBF(PI, PDI, PM1, PDM1, PM2, PDM2, NGL, ALPHA, BETA, ZG[I])
            JACOBF(PJ, PDJ, PM1, PDM1, PM2, PDM2, NGL, ALPHA, BETA, ZGL[J])
            FACI   = ALPHA*(ONE+ZG[I])-BETA*(ONE-ZG[I])
            FACJ   = ALPHA*(ONE+ZGL[J])-BETA*(ONE-ZGL[J])
            CONST  = EIGVAL*PJ+FACJ*PDJ
            D[I,J] = ((EIGVAL*PI+FACI*PDI)*(ZG[I]-ZGL[J])-
                     (ONE-ZG[I]^2)*PDI)/(CONST*(ZG[I]-ZGL[J])^2)
         end
         DT[J,I] = D[I,J]
      end
      end
      return nothing
      end
#
      function IGLM(I12, IT12, Z1, Z2, lz1, lz2, ND1, ND2)
#----------------------------------------------------------------------
#
#     Compute the one-dimensional interpolation operator (matrix) I12
#     ands its transpose IT12 for interpolating a variable from a
#     Gauss Legendre mesh (1) to a another mesh M (2).
#     Z1 : lz1 Gauss Legendre points.
#     Z2 : lz2 points on mesh M.
#
#--------------------------------------------------------------------
#      REAL I12[ND2,ND1],IT12[ND1,ND2],Z1[ND1],Z2[ND2]
      if lz1 == 1
         I12[1,1] = 1.
         IT12[1,1] = 1.
         return nothing
      end
      for I = 1:lz2
         ZI = Z2[I]
         for J = 1:lz1
            I12[I,J] = HGL(J, ZI, Z1, lz1)
            IT12[J,I] = I12[I,J]
      end
      end
      return nothing
      end
#
      function IGLLM(I12, IT12, Z1, Z2, lz1, lz2, ND1, ND2)
#----------------------------------------------------------------------
#
#     Compute the one-dimensional interpolation operator (matrix) I12
#     ands its transpose IT12 for interpolating a variable from a
#     Gauss-Lobatto Legendre mesh (1) to a another mesh M (2).
#     Z1 : lz1 Gauss-Lobatto Legendre points.
#     Z2 : lz2 points on mesh M.
#
#--------------------------------------------------------------------
#      REAL I12[ND2,ND1],IT12[ND1,ND2],Z1[ND1],Z2[ND2]
      if lz1 == 1
         I12[1,1] = 1.
         IT12[1,1] = 1.
         return nothing
      end
      for I = 1:lz2
         ZI = Z2[I]
         for J = 1:lz1
            I12[I,J] = HGLL(J, ZI, Z1, lz1)
            IT12[J,I] = I12[I,J]
      end
      end
      return nothing
      end
#
      function IGJM(I12, IT12, Z1, Z2, lz1, lz2, ND1, ND2, ALPHA, BETA)
#----------------------------------------------------------------------
#
#     Compute the one-dimensional interpolation operator (matrix) I12
#     ands its transpose IT12 for interpolating a variable from a
#     Gauss Jacobi mesh (1) to a another mesh M (2).
#     Z1 : lz1 Gauss Jacobi points.
#     Z2 : lz2 points on mesh M.
#     Single precision version.
#
#--------------------------------------------------------------------
#      REAL I12[ND2,ND1],IT12[ND1,ND2],Z1[ND1],Z2[ND2]
      if lz1 == 1
         I12[1,1] = 1.
         IT12[1,1] = 1.
         return nothing
      end
      for I = 1:lz2
         ZI = Z2[I]
         for J = 1:lz1
            I12[I,J] = HGJ(J, ZI, Z1, lz1, ALPHA, BETA)
            IT12[J,I] = I12[I,J]
      end
      end
      return nothing
      end
#
      function IGLJM(I12, IT12, Z1, Z2, lz1, lz2, ND1, ND2, ALPHA, BETA)
#----------------------------------------------------------------------
#
#     Compute the one-dimensional interpolation operator (matrix) I12
#     ands its transpose IT12 for interpolating a variable from a
#     Gauss-Lobatto Jacobi mesh (1) to a another mesh M (2).
#     Z1 : lz1 Gauss-Lobatto Jacobi points.
#     Z2 : lz2 points on mesh M.
#     Single precision version.
#
#--------------------------------------------------------------------
#      REAL I12[ND2,ND1],IT12[ND1,ND2],Z1[ND1],Z2[ND2]
      if lz1 == 1
         I12[1,1] = 1.
         IT12[1,1] = 1.
         return nothing
      end
      for I = 1:lz2
         ZI = Z2[I]
         for J = 1:lz1
            I12[I,J] = HGLJ(J, ZI, Z1, lz1, ALPHA, BETA)
            IT12[J,I] = I12[I,J]
      end
      end
      return nothing
      end
