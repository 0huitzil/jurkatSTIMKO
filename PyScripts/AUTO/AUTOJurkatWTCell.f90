!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Model2 :   Calcium exchange between the cytoplasm, ER, mitochondria and external medium via the PM
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PROGRAM model2
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)
      ! Variables 
      DOUBLE PRECISION c, ce, cm, h, Ct, s, i
      !Fluxes 
      DOUBLE PRECISION Jipr, JsercaC, JsercaM, Jdiff, cdiinfty
      DOUBLE PRECISION  Jh, Th, Jcrac
      ! Volume Parameters
      DOUBLE PRECISION delta, gammaE, gammaM
      ! IPR Parameters 
      DOUBLE PRECISION p, Kf, kbeta, Kp, Kc, Kh, Tmax, Kt
      DOUBLE PRECISION ma, A, B, ha, alpha, beta, P0 
      ! SERCA parameters 
      DOUBLE PRECISION VsC, VsM, KbarC, KbarM, Tg, Jpm, KsM, KsC
      ! diffusion parameters 
      DOUBLE PRECISION kdiff
      ! CDI parameters 
      DOUBLE PRECISION Ki, si, Ti
      ! CRAC parameters 
      DOUBLE PRECISION w1, w2, w12, w21, K1, K2, K12, K21, n, phi1, phi2, phi3, phi4, S1a, S2a, S12, S21, a0
      ! PM parameters 
      DOUBLE PRECISION Vplc, Vdeg, Kdeg, Tp, s1, Vsoce, Ts, Ke, Vpm, Kpm
      ! K1 parameters 
      DOUBLE PRECISION Vk, k1Rest, kn, nk, Tkmax, taun, kInfty, TauK
      DOUBLE PRECISION IPdeg, eps, PERIOD, TIME, num1, num2, num3
      !======================================
      ! Variables 
      c=U(1)
      ce=U(2)
      cm=U(3)
      h=U(4)
      p=U(5)
      s=U(6)
      i=U(7)
      ! K1=U(8)
      !======================================
      !Parameters 
      Vpm=PAR(1)
      Kpm=PAR(2)
      s1=PAR(3)
      Vsoce=PAR(4)
      Ts=PAR(5)
      Ke=PAR(6)
      a0=PAR(7)
      n=PAR(8)
      K1=PAR(9)
      K2=PAR(10)
      PERIOD=PAR(11)
      K12=PAR(12)
      K21=PAR(13)
      TIME=PAR(14)
      w1=PAR(15)
      w2=PAR(16)
      w12=PAR(17)
      w21=PAR(18)
      Vdeg=PAR(19)
      Kdeg=PAR(20)
      Vplc=PAR(21)
      Tp=PAR(22)
      gammaE=PAR(23)
      gammaM=PAR(24)
      delta=PAR(25)
      Ct=PAR(26)
      Kf=PAR(27)
      kbeta=PAR(28)
      Kp=PAR(29)
      Kc=PAR(30)
      Kh=PAR(31)
      Tmax=PAR(32)
      Kt=PAR(33)
      VsC=PAR(34)
      VsM=PAR(35)
      KbarC=PAR(36)
      KbarM=PAR(37)
      KsM=PAR(38)
      KsC=PAR(39)
      Tg=PAR(40)
      kdiff=PAR(41)
      Ki=PAR(42)
      si=PAR(43)
      Ti=PAR(44)
      Vk=PAR(45)
      k1Rest=PAR(46)
      kn=PAR(47)
      nk=PAR(48)
      Tkmax=PAR(49)
      taun=PAR(50)
      ! #======================================
      ! #IPR flux
      ma = c**4 / (Kc**4 + c**4)
      B = p**2 / (Kp**2 + p**2)
      A = 1 - p**2 / (Kp**2 + p**2)
      ha = Kh**4 / (Kh**4 + c**4)
      alpha = A*(1-ma*ha)
      beta = B*ma*h
      P0 = beta/(beta + kbeta*(beta + alpha))   
      Jipr=Kf*P0*(ce-c)
      ! #======================================
      ! #h fluxes
      Jh=(Kh**4 / (Kh**4 + c**4)) - h
      Th=Tmax*(Kt**4 / (Kt**4 + c**4))
      ! #======================================
      ! #SERCA flux
      JsercaC=(VsC*(c**2*(1-Tg) - KbarC*ce**2))/(c**2 + KsC**2)
      JsercaM=(VsM*(cm**2*(1-Tg) - KbarM*ce**2))/(cm**2 + KsM**2)
      ! #======================================
      ! #IP3 flux
      IPdeg=Vdeg*c**2/(c**2+Kdeg**2)
      ! #======================================
      ! #PMCA flux
      Jpm = (Vpm*c**2)/(c**2 + Kpm**2)
      ! #======================================
      ! #CRAC flux
      phi1=exp(n*(-K1 + ce))
      phi2=exp(n*(-K2 + ce))
      phi3=exp(n*(-K1 + K21 + ce))
      phi4=exp(K12*n)
      num1=(phi1 + 1)*(phi2 + 1)*(phi1*phi2 + phi1 + phi2 + 4*phi3 + 4*phi4 + 1)
      num2=2*(phi2*phi3 + phi2*phi4 + phi3 + phi4)
      num3= (phi1 + 1)/(2*(phi3 + phi4))
      S2a=sqrt(num1)/num2 - num3
      S1a=((phi2+1)/(phi1+1))*S2a
      S21=S1a*S2a*phi3
      S12=S1a*S2a*phi4
      Jcrac = w1*S1a + w2*S2a + w12*S12 + w21*S21
      ! #======================================
      ! #Linear diffusion
      Jdiff = kdiff*(cm-c)
      ! #======================================
      ! #CDI 
      cdiInfty = 1/(1+exp(si*(cm-Ki)))
      ! #======================================
      ! #Dynamic K1
      kInfty = (Vk*kn**8)/(kn**8 + ce**8) + k1Rest
      TauK=Tkmax*(ce**8 / (taun**8 + ce**8))
      ! #======================================
      !c
      F(1)=(Jipr - JsercaC + Jdiff) + delta*(-Jpm)
      !ce
      F(2)=gammaE*(JsercaC + JsercaM - Jipr)
      !cm
      F(3)=gammaM*(delta*s - JsercaM - Jdiff)
      !h
      F(4)=Jh/Th
      !p 
      F(5)=(Vplc - IPdeg*p)/Tp
      !s 
      F(6)=(Jcrac - s)/Ts
      !i
      ! F(7)=0
      !K1
      ! F(8)=(kInfty-K1)/Tauk 
END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      ! Variables 
      DOUBLE PRECISION c, ce, cm, h, Ct, s, i
      !Fluxes 
      DOUBLE PRECISION Jipr, JsercaC, JsercaM, Jdiff, cdiinfty
      DOUBLE PRECISION  Jh, Th, Jcrac
      ! Volume Parameters
      DOUBLE PRECISION delta, gammaE, gammaM
      ! IPR Parameters 
      DOUBLE PRECISION p, Kf, kbeta, Kp, Kc, Kh, Tmax, Kt
      DOUBLE PRECISION ma, A, B, ha, alpha, beta, P0 
      ! SERCA parameters 
      DOUBLE PRECISION VsC, VsM, KbarC, KbarM, Tg, Jpm, KsM, KsC
      ! diffusion parameters 
      DOUBLE PRECISION kdiff
      ! CDI parameters 
      DOUBLE PRECISION Ki, si, Ti
      ! CRAC parameters 
      DOUBLE PRECISION w1, w2, w12, w21, K1, K2, K12, K21, n, phi1, phi2, phi3, phi4, S1a, S2a, S12, S21, a0
      ! PM parameters 
      DOUBLE PRECISION Vplc, Vdeg, Kdeg, Tp, s1, Vsoce, Ts, Ke, Vpm, Kpm
      ! K1 parameters 
      DOUBLE PRECISION Vk, k1Rest, kn, nk, Tkmax, taun, kInfty, TauK
      DOUBLE PRECISION IPdeg, eps, PERIOD, TIME, num1, num2, num3
      !======================================
      ! # PMCA 
      Vpm=1.0
      Kpm=0.2
      ! #======================================
      ! # CRAC 
      s1=0.15
      Vsoce=3
      Ts=4
      Ke=800
      n=0.15
      a0=0
      K1=200
      K2=350
      K12=400
      K21=470
      w1=1
      w2=0.7
      w12=2
      w21=2
      ! #======================================
      ! # IP3
      ! #Vplc=0.01 for burst, 0.005 for spike
      Vdeg=6
      Kdeg=0.5
      Vplc=0
      Tp=2.9
      ! #======================================
      ! # General 
      gammaE=5.5
      gammaM=50
      delta=3
      Ct=200
      ! #======================================
      ! # IPR 
      Kf=1.6
      kbeta=0.4
      Kp=10
      Kc=0.16
      ! #======================================
      ! # H 
      Kh=0.168
      Tmax=10
      Kt=0.095
      ! #======================================
      ! # SERCA 
      VsC=1
      VsM=0.2
      KbarC=2.2e-08
      KbarM=2.2e-08
      KsM=0.2
      KsC=0.2
      Tg=0
      ! #======================================
      ! # Linear diffusion 
      kdiff=0.004
      ! #======================================
      ! # CDI
      Ki=800
      si=0.2
      Ti=1
      ! #======================================
      ! # Dynamic K1
      Vk=670
      k1Rest=200
      kn=100
      nk=8
      Tkmax=800
      taun=60
      !======================================
      PAR(1)=Vpm
      PAR(2)=Kpm
      PAR(3)=s1
      PAR(4)=Vsoce
      PAR(5)=Ts
      PAR(6)=Ke
      PAR(7)=a0
      PAR(8)=n
      PAR(9)=K1
      PAR(10)=K2
      PAR(11)=0.0
      PAR(12)=K12
      PAR(13)=K21
      PAR(14)=0.0
      PAR(15)=w1
      PAR(16)=w2
      PAR(17)=w12
      PAR(18)=w21
      PAR(19)=Vdeg
      PAR(20)=Kdeg
      PAR(21)=Vplc
      PAR(22)=Tp
      PAR(23)=gammaE
      PAR(24)=gammaM
      PAR(25)=delta
      PAR(26)=Ct
      PAR(27)=Kf
      PAR(28)=kbeta
      PAR(29)=Kp
      PAR(30)=Kc
      PAR(31)=Kh
      PAR(32)=Tmax
      PAR(33)=Kt
      PAR(34)=VsC
      PAR(35)=VsM
      PAR(36)=KbarC
      PAR(37)=KbarM
      PAR(38)=KsM
      PAR(39)=KsC
      PAR(40)=Tg
      PAR(41)=kdiff
      PAR(42)=Ki
      PAR(43)=si
      PAR(44)=Ti
      PAR(45)=Vk
      PAR(46)=k1Rest
      PAR(47)=kn
      PAR(48)=nk
      PAR(49)=Tkmax
      PAR(50)=taun 
      !======================================
      U(1)=7.86161E-02
      U(2)=8.37105E+02
      U(3)=5.04548E+01
      U(4)=9.54242E-01
      U(5)=0
      U(6)=1.33833E-01
      ! U(7)=1
      ! U(8)=2.00762E+02
END SUBROUTINE STPNT

SUBROUTINE BCND
END SUBROUTINE BCND

SUBROUTINE ICND
END SUBROUTINE ICND

SUBROUTINE FOPT
END SUBROUTINE FOPT

SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      ! Variables 
      DOUBLE PRECISION c, ce, cm, h, Ct, s
      !Fluxes 
      DOUBLE PRECISION Jipr, JsercaC, JsercaM, Jdiff
      DOUBLE PRECISION  Jh, Th
      ! Volume Parameters
      DOUBLE PRECISION delta, gammaE, gammaM
      ! IPR Parameters 
      DOUBLE PRECISION p, Kf, kbeta, Kp, Kc, Kh, Tmax, Kt
      DOUBLE PRECISION ma, A, B, ha, alpha, beta, P0 
      ! SERCA parameters 
      DOUBLE PRECISION Vs, Kbar, Tg, Jpm, VsM
      ! diffusion parameters 
      DOUBLE PRECISION kdiff
      ! PM parameters 
      DOUBLE PRECISION Vplc, Vdeg, Kdeg, Tp, s1, Vsoce, Ts, Ke, Vpm, Kpm
      DOUBLE PRECISION IPdeg, L, Jsoce, eps, PERIOD, TIME
      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      !======================================
      ! Frequency 
      PAR(99)=1/PAR(11)
      !======================================
      ! Variables 
      c=U(1)
      ce=U(2)
      cm=U(3)
      h=U(4)
      p=U(5)
      s=U(6)
      
      Vdeg=PAR(15)
      Kdeg=PAR(16)
      Vplc=PAR(17)
      IPdeg=Vdeg*c**2/(c**2+Kdeg**2)
      L=Vplc - IPdeg*p
      PAR(98)=GETP('MIN',8,U)
      PAR(97)=GETP('MIN',7,U)
      PAR(96)=GETP('MIN',6,U)
END SUBROUTINE PVLS

! END PROGRAM model2