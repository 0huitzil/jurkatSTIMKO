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
      ! i=U(7)
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
      ! p=PAR(26)
      ! ce = gamma*(Ct-c)
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
      Jcrac = w2/(1+exp(n*(ce-K2)))
      ! #======================================
      ! #Linear diffusion
      Jdiff = kdiff*(cm-c)
      ! #======================================
      ! #CDI 
      ! cdiInfty = 1/(1+exp(si*(cm-Ki)))
      ! #======================================
      ! #Dynamic K1
      ! kInfty = (Vk*kn**8)/(kn**8 + ce**8) + k1Rest
      ! TauK=Tkmax*(ce**8 / (taun**8 + ce**8))
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
      ! F(7)=(cdiInfty - i)/Ti
      !K1
      ! F(8)=(kInfty-K1)/Tauk 
      IF(IJAC.EQ.0)RETURN

      DFDU(1,1)= -4*Kf*c**7*h*p**2*(-c + ce)/((Kc**4 + c**4)**2*(Kp**2 &  
      + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2))  & 
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) &  
      + Kf*c**4*h*p**2*(-c + ce)*(4*c**7*h*p**2/((Kc**4 + c**4)**2*(Kp**2 + p**2)) & 
      - 4*c**3*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      - kbeta*(-4*c**7*h*p**2/((Kc**4 + c**4)**2*(Kp**2 + p**2)) &  
      + 4*c**3*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + (-p**2/(Kp**2 + p**2) + 1)*(4*Kh**4*c**7/((Kc**4 + c**4)*(Kh**4 + c**4)**2) &  
      + 4*Kh**4*c**7/((Kc**4 + c**4)**2*(Kh**4 + c**4)) & 
      - 4*Kh**4*c**3/((Kc**4 + c**4)*(Kh**4 + c**4)))))/((Kc**4 + c**4)*(Kp**2 &  
      + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))**2) & 
      - Kf*c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) & 
      + 4*Kf*c**3*h*p**2*(-c + ce)/((Kc**4 + c**4)*(Kp**2 &  
      + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) & 
      + 2*Vpm*c**3*delta/(Kpm**2 + c**2)**2 - 2*Vpm*c*delta/(Kpm**2 + c**2) &  
      - 2*VsC*c*(1 - Tg)/(KsC**2 + c**2) & 
      + 2*VsC*c*(-KbarC*ce**2 + c**2*(1 - Tg))/(KsC**2 + c**2)**2 - kdiff
      DFDU(1,2)= 2*KbarC*VsC*ce/(KsC**2 + c**2) &  
      + Kf*c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1))))
      DFDU(1,3)= kdiff
      DFDU(1,4)= Kf*c**4*h*p**2*(-c + ce)*(-c**4*kbeta*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      - c**4*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)))/((Kc**4 & 
      + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))**2) & 
      + Kf*c**4*p**2*(-c + ce)/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2))  &
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1))))
      DFDU(1,5)= -2*Kf*c**4*h*p**3*(-c + ce)/((Kc**4 + c**4)*(Kp**2 + p**2)**2*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2))&  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) & 
      + Kf*c**4*h*p**2*(-c + ce)*(2*c**4*h*p**3/((Kc**4 + c**4)*(Kp**2 + p**2)**2) & 
      - 2*c**4*h*p/((Kc**4 + c**4)*(Kp**2 + p**2)) & 
      - kbeta*(-2*c**4*h*p**3/((Kc**4 + c**4)*(Kp**2 + p**2)**2) & 
      + 2*c**4*h*p/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (2*p**3/(Kp**2 + p**2)**2 - 2*p/(Kp**2 + p**2))*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) &  
      + 1)))/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))**2) &  
      + 2*Kf*c**4*h*p*(-c + ce)/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1))))
      DFDU(1,6)= 0
      DFDU(2,1)= gammaE*(4*Kf*c**7*h*p**2*(-c + ce)/((Kc**4 &  
      + c**4)**2*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) &  
      - Kf*c**4*h*p**2*(-c + ce)*(4*c**7*h*p**2/((Kc**4 + c**4)**2*(Kp**2 + p**2)) &  
      - 4*c**3*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      - kbeta*(-4*c**7*h*p**2/((Kc**4 + c**4)**2*(Kp**2 + p**2))&  
      + 4*c**3*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(4*Kh**4*c**7/((Kc**4 + c**4)*(Kh**4 + c**4)**2) &  
      + 4*Kh**4*c**7/((Kc**4 + c**4)**2*(Kh**4 + c**4)) &  
      - 4*Kh**4*c**3/((Kc**4 + c**4)*(Kh**4 + c**4)))))/((Kc**4 + c**4)*(Kp**2 &  
      + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))**2) &  
      + Kf*c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 &  
      + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) &  
      - 4*Kf*c**3*h*p**2*(-c + ce)/((Kc**4 + c**4)*(Kp**2 &  
      + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) &  
      + 2*VsC*c*(1 - Tg)/(KsC**2 + c**2) &  
      - 2*VsC*c*(-KbarC*ce**2 + c**2*(1 - Tg))/(KsC**2 + c**2)**2)
      DFDU(2,2)= gammaE*(-2*KbarC*VsC*ce/(KsC**2 + c**2) &  
      - 2*KbarM*VsM*ce/(KsM**2 + cm**2) &  
      - Kf*c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 &  
      + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))))
      DFDU(2,3)= gammaE*(2*VsM*cm*(1 - Tg)/(KsM**2 + cm**2) &  
      - 2*VsM*cm*(-KbarM*ce**2 + cm**2*(1 - Tg))/(KsM**2 + cm**2)**2)
      DFDU(2,4)= gammaE*(-Kf*c**4*h*p**2*(-c + ce)*(-c**4*kbeta*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      - c**4*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)))/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))**2) &  
      - Kf*c**4*p**2*(-c + ce)/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))))
      DFDU(2,5)= gammaE*(2*Kf*c**4*h*p**3*(-c + ce)/((Kc**4 + c**4)*(Kp**2 + p**2)**2*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))) &  
      - Kf*c**4*h*p**2*(-c + ce)*(2*c**4*h*p**3/((Kc**4 + c**4)*(Kp**2 + p**2)**2) - 2*c**4*h*p/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      - kbeta*(-2*c**4*h*p**3/((Kc**4 + c**4)*(Kp**2 + p**2)**2) + 2*c**4*h*p/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (2*p**3/(Kp**2 + p**2)**2 - 2*p/(Kp**2 + p**2))*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) &  
      + 1)))/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))**2) &  
      - 2*Kf*c**4*h*p*(-c + ce)/((Kc**4 + c**4)*(Kp**2 + p**2)*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + kbeta*(c**4*h*p**2/((Kc**4 + c**4)*(Kp**2 + p**2)) &  
      + (-p**2/(Kp**2 + p**2) + 1)*(-Kh**4*c**4/((Kc**4 + c**4)*(Kh**4 + c**4)) + 1)))))
      DFDU(2,6)= 0
      DFDU(3,1)= gammaM*kdiff
      DFDU(3,2)= 2*KbarM*VsM*ce*gammaM/(KsM**2 + cm**2)
      DFDU(3,3)= gammaM*(-2*VsM*cm*(1 - Tg)/(KsM**2 + cm**2) &  
      + 2*VsM*cm*(-KbarM*ce**2 + cm**2*(1 - Tg))/(KsM**2 + cm**2)**2 - kdiff)
      DFDU(3,4)= 0
      DFDU(3,5)= 0
      DFDU(3,6)= delta*gammaM
      DFDU(4,1)= -4*Kh**4*c**3*(Kt**4 + c**4)/(Kt**4*Tmax*(Kh**4 + c**4)**2) &  
      + 4*c**3*(Kh**4/(Kh**4 + c**4) - h)/(Kt**4*Tmax)
      DFDU(4,2)= 0
      DFDU(4,3)= 0
      DFDU(4,4)= -(Kt**4 + c**4)/(Kt**4*Tmax)
      DFDU(4,5)= 0
      DFDU(4,6)= 0
      DFDU(5,1)= (2*Vdeg*c**3*p/(Kdeg**2 + c**2)**2 - 2*Vdeg*c*p/(Kdeg**2 + c**2))/Tp
      DFDU(5,2)= 0
      DFDU(5,3)= 0
      DFDU(5,4)= 0
      DFDU(5,5)= -Vdeg*c**2/(Tp*(Kdeg**2 + c**2))
      DFDU(5,6)= 0
      DFDU(6,1)= 0
      DFDU(6,2)= -n*w2*exp(n*(-K2 + ce))/(Ts*(exp(n*(-K2 + ce)) + 1)**2)
      DFDU(6,3)= 0
      DFDU(6,4)= 0
      DFDU(6,5)= 0
      DFDU(6,6)= -1/Ts

      IF(IJAC.EQ.1)RETURN
      ! Only parameter derivatives for the Vplc parameter are provided
      DFDP(1,21)= 0
      DFDP(2,21)= 0
      DFDP(3,21)= 0
      DFDP(4,21)= 0
      DFDP(5,21)= 1/Tp
      DFDP(6,21)= 0
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
      U(1)=2.78105E-02
      U(2)=3.73872E+02
      U(3)=1.41474E-01
      U(4)=9.99250E-01
      U(5)=0
      U(6)=1.89689E-02
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
      double precision l1r, l1i, l2r, l2i, l3r ,l3i ,l4r ,l4i ,l5r
      DOUBLE PRECISION l5i ,l6r ,l6i ,l7r ,l7i
      DOUBLE PRECISION l1, l2, l3, l4, l5, l6, l7
      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      ! The first argument of GETP may be one of the following:
      !        'NRM' (L2-norm),     'MAX' (maximum),
      !        'INT' (integral),    'BV0 (left boundary value),
      !        'MIN' (minimum),     'BV1' (right boundary value).
      !        'MNT' (t value for minimum)
      !        'MXT' (t value for maximum)
      !        'NDIM', 'NDX' (effective (active) number of dimensions)
      !        'NTST' (NTST from constant file)
      !        'NCOL' (NCOL from constant file)
      !        'NBC'  (active NBC)
      !        'NINT' (active NINT)
      !        'DTM'  (delta t for all t values, I=1...NTST)
      !        'WINT' (integration weights used for interpolation, I=0...NCOL)
      !
      ! Also available are
      !   'STP' (Pseudo-arclength step size used).
      !   'FLD' (`Fold function', which vanishes at folds).
      !   'BIF' (`Bifurcation function', which vanishes at singular points).
      !   'HBF' (`Hopf function'; which vanishes at Hopf points).
      !   'SPB' ( Function which vanishes at secondary periodic bifurcations).
      !   'EIG' ( Eigenvalues/multipliers, I=1...2*NDIM, alternates real/imag parts).
      !   'STA' ( Number of stable eigenvalues/multipliers).
      ! Frequency 
      PAR(99)=1/PAR(11)
      !======================================
      ! Variables 
      ! c=U(1)
      ! ce=U(2)
      ! cm=U(3)
      ! h=U(4)
      ! p=U(5)
      ! s=U(6)
      ! !======================================
      ! ! Multipliers/eigenvalues
      ! l1r = GETP('EIG',1,U)
      ! l1i = GETP('EIG',2,U)
      ! l2r = GETP('EIG',3,U)
      ! l2i = GETP('EIG',4,U)
      ! l3r = GETP('EIG',5,U)
      ! l3i = GETP('EIG',6,U)
      ! l4r = GETP('EIG',7,U)
      ! l4i = GETP('EIG',8,U)
      ! l5r = GETP('EIG',9,U)
      ! l5i = GETP('EIG',10,U)
      ! l6r = GETP('EIG',11,U)
      ! l6i = GETP('EIG',12,U)
      ! l7r = GETP('EIG',13,U)
      ! l7i = GETP('EIG',14,U)
      ! ! Modulus
      ! ! l1 = sqrt(l1r**2 + l1i**2)
      ! l2 = sqrt(l2r**2 + l2i**2)
      ! l3 = sqrt(l3r**2 + l3i**2)
      ! l4 = sqrt(l4r**2 + l4i**2)
      ! l5 = sqrt(l5r**2 + l5i**2)
      ! l6 = sqrt(l6r**2 + l6i**2)
      ! l7 = sqrt(l7r**2 + l7i**2)
      ! PAR(80)=l2
      ! PAR(81)=l3
      ! PAR(82)=l4
      ! PAR(83)=l5
      ! PAR(84)=l6
      ! PAR(85)=l7
      ! PAR(98)=GETP('MIN',8,U)
      ! PAR(97)=GETP('MIN',7,U)
      ! PAR(96)=GETP('MIN',6,U)
END SUBROUTINE PVLS

! END PROGRAM model2