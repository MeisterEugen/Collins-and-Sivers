!Fortran module for the calculation of spin-dependent quantities.
MODULE ROUTINES
!Module with routines to handle four vectors.
use NewVariables

!kind for real variables
integer, parameter :: rk=selected_real_kind(8)

!Parameters for 3P0 mechanism: complex mass.
complex(kind=rk) :: mu = (0.42,0.76)

!Parameters for vector meson production.
!Ratio of couplings to vector mesons with L and T polarization
!w.r.t the string axis.
real(kind=rk)    :: GLGT = 1 
!Oblique polarization of vector mesons.
real(kind=rk)    :: THETALT = 0*acos(-1.0_rk)/2
!Form factor (given by Breit-Wigner functions) in the 3-body decays
!of omega and phi mesons.
logical :: FormFact = .false.

!Some useful numbers.
real(kind=rk)       :: pi = acos(-1.)
complex(kind=rk)    :: im = (0._rk,1._rk),cu=(1._rk,0._rk)

!Definition of Pauli matrices.
REAL,DIMENSION(2,2)     :: SIG1=transpose(reshape((/0,1,1,0/),shape(SIG1)))
REAL,DIMENSION(2,2)     :: SIG3=transpose(reshape((/1,0,0,-1/),shape(SIG3)))
REAL,DIMENSION(2,2)     :: IDEN=transpose(reshape((/1,0,0,1/),shape(IDEN)))
COMPLEX,DIMENSION(2,2)  :: SIG2=transpose(reshape((/(0._rk,0._rk),(0._rk,-1._rk),(0._rk,1._rk),(0._rk,0._rk)/),shape(SIG2)))

!Masses of hadrons in GeV, needed for the decays.
!Pseudoscalar mesons.
real(kind=rk) :: Mpi    = 0.139570
real(kind=rk) :: Mpiz   = 0.134977
real(kind=rk) :: MK     = 0.493667
real(kind=rk) :: MKz    = 0.497648
real(kind=rk) :: Meta   = 0.547853
real(kind=rk) :: Metap  = 0.95766
real(kind=rk) :: MKzL   = 0.497648
real(kind=rk) :: MKzS   = 0.497648
!Vector mesons.
real(kind=rk) :: Mrho       = 0.77511
real(kind=rk) :: Mrhoz      = 0.77526
real(kind=rk) :: Momega     = 0.78265
real(kind=rk) :: Mphi       = 1.019455
real(kind=rk) :: MKstar     = 0.89166
real(kind=rk) :: MKstarz    = 0.89581
!Baryons.
real(kind=rk) :: MN         = 0.939565
real(kind=rk) :: MP         = 0.938272
real(kind=rk) :: Mdeuterium = 1.875612

!Decay widths of vector mesons.
real(kind=rk) :: Gammarhozpipi      = 0.1462
real(kind=rk) :: Gammarhoppipi      = 0.1491
real(kind=rk) :: Gammarhompipi      = 0.1491
real(kind=rk) :: GammaKstarpKpi     = 0.0508
real(kind=rk) :: GammaKstarmKpi     = 0.0508
real(kind=rk) :: GammaKstarzKzpi    = 0.0474
real(kind=rk) :: GammaKstarzbKzbpi  = 0.0474
real(kind=rk) :: GammaPhi           = 0.00426
real(kind=rk) :: GammaOmega         = 0.00849

!Flavor codes of pseudoscalar mesons.
integer :: KFpip    = 211
integer :: KFpim    = -211
integer :: KFpiz    = 111
integer :: KFKp     = 321
integer :: KFKm     = -321
integer :: KFKz     = 311
integer :: KFKzb    = -311
integer :: KFeta    = 221
integer :: KFetap   = 331
integer :: KFKzL    = 130
integer :: KFKzS    = 310

!Flavor codes of vector mesons.
integer :: KFrhop    = 213
integer :: KFrhom    = -213
integer :: KFrhoz    = 113
integer :: KFKstarp  = 323
integer :: KFKstarm  = -323
integer :: KFKstarz  = 313
integer :: KFKstarzb = -313
integer :: KFomega   = 223
integer :: KFphi     = 333

!Flavor codes for quarks.
integer :: KFLu = 2
integer :: KFLd = 1
integer :: KFLs = 3

!Flavour codes for vector bosons.
integer :: KFphoton = 22

!Masses of vector bosons.
real(kind=rk) :: Mphoton = 0.0

CONTAINS

!Calculates the polarization vector of q' in q -> PS + q'.
SUBROUTINE OUTPOL(KTPRIMX,KTPRIMY,SX,SY,SZ,SOUTX,SOUTY,SOUTZ)
	IMPLICIT NONE
	REAL(KIND=RK), INTENT(IN) :: KTPRIMX,KTPRIMY
	REAL(KIND=RK), INTENT(IN) :: SX, SY, SZ
	REAL(KIND=RK), INTENT(OUT) :: SOUTX, SOUTY, SOUTZ
	REAL(KIND=RK), DIMENSION(2) :: KTPRIM
	REAL(KIND=RK), DIMENSION(3) :: S
	COMPLEX(kind=rk), DIMENSION(2,2)::RHO,V1,V2,RHOV1,V2RHOV1,MATX,MATY,MATZ
	
    !Transverse momentum of q'.
	KTPRIM(1)=KTPRIMX
	KTPRIM(2)=KTPRIMY
	!Polarization vector of q.
	S(1)=SX
	S(2)=SY
	S(3)=SZ
	!Density matrix of q.
	RHO = 0.5_rk*(IDEN+S(1)*SIG1+S(2)*SIG2+S(3)*SIG3)
	!Vertex (3P0 propagator x sigma_z) q-PS-q' and its complex conjugate.
	V1 = CONJG(MU)*SIG3 - KTPRIM(1)*SIG1 - KTPRIM(2)*SIG2
	V2 = MU*SIG3 - KTPRIM(1)*SIG1 - KTPRIM(2)*SIG2
			
	!Calculation of the density matrix of q'.
	RHOV1 = MATMUL(RHO,V1)
	V2RHOV1 = MATMUL(V2,RHOV1)
	RHO = V2RHOV1/(V2RHOV1(1,1) + V2RHOV1(2,2))
	MATX=MATMUL(SIG1,RHO)
	MATY=MATMUL(SIG2,RHO)
	MATZ=MATMUL(SIG3,RHO)
	! Polarization vector of q'.
	SOUTX=(MATX(1,1)+MATX(2,2))
	SOUTY=(MATY(1,1)+MATY(2,2))
	SOUTZ=(MATZ(1,1)+MATZ(2,2))
	
END SUBROUTINE

!Calculates the polarization vector of q' in q -> VM + q', assuming no
!information is coming from the VM decay process (no decay matrix).
SUBROUTINE OUTPOL_VM_noDecayMatrix(KTPRIMX,KTPRIMY,SX,SY,SZ,SOUTX,SOUTY,SOUTZ)
    IMPLICIT NONE
    REAL(KIND=RK), INTENT(IN) :: KTPRIMX,KTPRIMY
    REAL(KIND=RK), INTENT(IN) :: SX, SY, SZ
    REAL(KIND=RK), INTENT(OUT) :: SOUTX, SOUTY, SOUTZ
    REAL(KIND=RK), DIMENSION(2) :: KTPRIM
    REAL(KIND=RK), DIMENSION(3) :: S
    COMPLEX(kind=rk), DIMENSION(2,2)::RHO,RHO0,DELTA,DELTADAG,MATX,MATY,MATZ
    COMPLEX(kind=rk), DIMENSION(2,2)::RHODELTADAG,DELTARHODELTADAG
    REAL(KIND=RK)  :: DDT, DDL
    COMPLEX(kind=rk)::trace
    
    !Transverse momentum of q'.
    KTPRIM(1)=KTPRIMX
    KTPRIM(2)=KTPRIMY
    !Polarization vector of q.
    S(1)=SX
    S(2)=SY
    S(3)=SZ
    !3P0 propagator.
    DELTA = MU*IDEN + MATMUL(SIG3,SIG1)*KTPRIM(1) + MATMUL(SIG3,SIG2)*KTPRIM(2)
    DELTADAG = CONJG(MU)*IDEN - MATMUL(SIG3,SIG1)*KTPRIM(1) - MATMUL(SIG3,SIG2)*KTPRIM(2)

    !Spin density matrix of q after multiplying by the VM coupling.
    DDT = (GLGT**2) / (2.0+GLGT**2)
    DDL = (GLGT**2-2.0) / (2.0+GLGT**2)
    RHO=0.5_rk*(IDEN + DDT*S(1)*SIG1 + DDT*S(2)*SIG2 + DDL*S(3)*SIG3)
    !Calculation of the spin density matrix of q'.
    RHO=MATMUL(RHO,DELTADAG)
    RHO=MATMUL(DELTA,RHO)
    RHO=RHO/(RHO(1,1)+RHO(2,2))

    !Polarization vector of q'.
    MATX = MATMUL(SIG1,RHO)
    MATY = MATMUL(SIG2,RHO)
    MATZ = MATMUL(SIG3,RHO)
    SOUTX=(MATX(1,1)+MATX(2,2))
    SOUTY=(MATY(1,1)+MATY(2,2))
    SOUTZ=(MATZ(1,1)+MATZ(2,2))
    
END SUBROUTINE

! Calculates the polarization vector of q' in q -> VM + q', taking into account
! the decay matrix coming from the decay of the VM.
SUBROUTINE OUTPOL_VM(KTPRIMX,KTPRIMY,SX,SY,SZ,SOUTX,SOUTY,SOUTZ,&
                     DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ)
    IMPLICIT NONE
    REAL(KIND=RK), INTENT(IN) :: KTPRIMX,KTPRIMY
    REAL(KIND=RK), INTENT(IN) :: SX, SY, SZ
    REAL(KIND=RK), INTENT(IN) :: DXX, DXY, DXZ
    REAL(KIND=RK), INTENT(IN) :: DYX, DYY, DYZ
    REAL(KIND=RK), INTENT(IN) :: DZX, DZY, DZZ
    REAL(KIND=RK), INTENT(OUT) :: SOUTX, SOUTY, SOUTZ
    REAL(KIND=RK), DIMENSION(2) :: KTPRIM
    REAL(KIND=RK), DIMENSION(3) :: S
    COMPLEX(kind=rk), DIMENSION(2,2)::RHO,RHO0,DELTA,DELTADAG,MATX,MATY,MATZ
    !Matrices for the calculation of the VM spin density matrix.
    COMPLEX(kind=rk), DIMENSION(2,2) :: XRHOX,XRHOY,XRHOZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: YRHOX,YRHOY,YRHOZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: ZRHOX,ZRHOY,ZRHOZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: MRHOXX,MRHOXY,MRHOXZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: MRHOYX,MRHOYY,MRHOYZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: MRHOZX,MRHOZY,MRHOZZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: RHOX_QP,RHOY_QP,RHOZ_QP
    COMPLEX(kind=rk)				 :: TRACE_QP
    
    !Transverse momentum of q'.
    KTPRIM(1)=KTPRIMX
    KTPRIM(2)=KTPRIMY

    !Polarization vector of q.
    S(1)=SX
    S(2)=SY
    S(3)=SZ

    !3P0 propagator.
    DELTA = MU*IDEN + MATMUL(SIG3,SIG1)*KTPRIM(1) + MATMUL(SIG3,SIG2)*KTPRIM(2)
    DELTADAG = CONJG(MU)*IDEN - MATMUL(SIG3,SIG1)*KTPRIM(1) - MATMUL(SIG3,SIG2)*KTPRIM(2)

    !Density matrix of q.
    RHO=0.5_rk*(IDEN + S(1)*SIG1 + S(2)*SIG2 + S(3)*SIG3)

    !Matrix multiplications.
    ! sigma_i x rho x sigma_j
    XRHOX = MATMUL(SIG1, MATMUL(RHO, SIG1))
    XRHOY = MATMUL(SIG1, MATMUL(RHO, SIG2))
    XRHOZ = MATMUL(SIG1, RHO)

    YRHOY = MATMUL(SIG2, MATMUL(RHO, SIG2))
    YRHOX = MATMUL(SIG2, MATMUL(RHO, SIG1))
    YRHOZ = MATMUL(SIG2, RHO)
    
    ZRHOZ = RHO
    ZRHOX = MATMUL(RHO, SIG1)
    ZRHOY = MATMUL(RHO, SIG2)

    !3P0_prpagator x [sigma_i x rho x sigma_j] x 3P0_propagator_dagger
    MRHOXX = MATMUL(DELTA,MATMUL(YRHOY,DELTADAG))
    MRHOXY = -MATMUL(DELTA,MATMUL(YRHOX,DELTADAG))
    MRHOXZ = -MATMUL(DELTA,MATMUL(YRHOZ,DELTADAG))*(GLGT)*(SIN(THETALT)+IM*COS(THETALT))

    MRHOYY = MATMUL(DELTA,MATMUL(XRHOX,DELTADAG))
    MRHOYX = -MATMUL(DELTA,MATMUL(XRHOY,DELTADAG))
    MRHOYZ = MATMUL(DELTA,MATMUL(XRHOZ,DELTADAG))*(GLGT)*(SIN(THETALT)+IM*COS(THETALT))
  
    MRHOZZ = MATMUL(DELTA,MATMUL(ZRHOZ,DELTADAG))*(GLGT**2)
    MRHOZX = MATMUL(DELTA,MATMUL(ZRHOY,DELTADAG))*(GLGT)*(-SIN(THETALT)+IM*COS(THETALT))
    MRHOZY = MATMUL(DELTA,MATMUL(ZRHOX,DELTADAG))*(GLGT)*(SIN(THETALT)-IM*COS(THETALT))

    !Calculation of the density matrix of q'.
	RHOX_QP = DXX * MRHOXX + DXY * MRHOXY + DXZ * MRHOXZ
	RHOY_QP = DYX * MRHOYX + DYY * MRHOYY + DYZ * MRHOYZ
	RHOZ_QP = DZX * MRHOZX + DZY * MRHOZY + DZZ * MRHOZZ
    TRACE_QP=RHOX_QP(1,1)+RHOX_QP(2,2)+RHOY_QP(1,1)+RHOY_QP(2,2)+RHOZ_QP(1,1)+RHOZ_QP(2,2)
	RHO=(RHOX_QP+RHOY_QP+RHOZ_QP)/TRACE_QP

    !Polarization vector of q'.
    MATX = MATMUL(SIG1,RHO)
    MATY = MATMUL(SIG2,RHO)
    MATZ = MATMUL(SIG3,RHO)
    SOUTX=(MATX(1,1)+MATX(2,2))
    SOUTY=(MATY(1,1)+MATY(2,2))
    SOUTZ=(MATZ(1,1)+MATZ(2,2))
    
END SUBROUTINE

!Calculates the polarization vector of q' in q -> PS/VM + q',
!without decay matrix in the VM case.
subroutine PropagateSpin(KTPRIMX,KTPRIMY,SX,SY,SZ,hadron,SOUTX,SOUTY,SOUTZ)
    implicit none
    real(kind=rk), intent(in) :: KTPRIMX,KTPRIMY
    real(kind=rk), intent(in) :: SX, SY, SZ
    integer, intent(in)       :: hadron
    real(kind=rk), intent(out) :: SOUTX, SOUTY, SOUTZ
    real(kind=rk), dimension(2) :: KTPRIM

    if(hadron.eq.0) then
        call OUTPOL(KTPRIMX,KTPRIMY,SX,SY,SZ,SOUTX,SOUTY,SOUTZ)
    else if(hadron.eq.1) then
        call OUTPOL_VM_noDecayMatrix(KTPRIMX,KTPRIMY,SX,SY,SZ,SOUTX,SOUTY,SOUTZ)
    else
        SOUTX = 0.0_rk
        SOUTY = 0.0_rk
        SOUTZ = 0.0_rk
    endif

end subroutine PropagateSpin

!Calculation of the 3P0 weight for PS or VM.
SUBROUTINE ACCEPTPYTHIA(KTPX,KTPY,SX,SY,FLAG,isVM,weight)
    IMPLICIT NONE
    REAL(KIND=RK), INTENT(IN) ::KTPX,KTPY,SX,SY
    INTEGER, INTENT(IN) :: isVM
    INTEGER,INTENT(OUT) :: FLAG
    REAL(KIND=RK),INTENT(OUT) :: weight
    REAL(KIND=RK) ::AP,W,RND, Coupling

	FLAG=1
    if(isVM .eq. 1) then
      Coupling = (GLGT**2)/(2.0+(GLGT**2))
    else
      Coupling = -1.0_rk
    endif

    AP=2.*AIMAG(MU)*Coupling/(REAL(MU)**2+AIMAG(MU)**2+KTPX**2+KTPY**2)
	W=(1.+AP*(KTPX*SY-KTPY*SX))!/2.0
	weight=W
	CALL RANDOM_NUMBER(RND)
	IF(RND<W) FLAG=0
if(weight/2>1.0 .or. weight/2<0) print*, weight, AP * (KTPX*SY-KTPY*SX), sqrt(SX**2+SY**2), isVM
	
END SUBROUTINE

!Set a default seed for the randum number generator.
subroutine setSeed()
    implicit none
    integer,dimension(33) :: seed=(/1769108563, 1414555502, 1752065338,&
                                  1701015137, 1684625252, 1025534068,&
                                  3026976, 123431, 389276, 273, 222991,&
                                  1769108563, 1414555502, 1752065338,&
                                  1701015137, 1684625252, 1025534068,&
                                  3026976, 123431, 389276, 273, 222991,&
                                  1769108563, 1414555502, 1752065338,&
                                  1701015137, 1684625252, 1025534068,&
                                  3026976, 123431, 389276, 273, 222991 /)
    call random_seed(put = seed)
end subroutine setSeed

!Setup for complex mass.
SUBROUTINE SETMU(RE,IM)
    IMPLICIT NONE
    REAL(KIND=RK),INTENT(IN)	:: RE,IM
    COMPLEX(KIND=RK)			:: NEW_MU

	NEW_MU = CMPLX(RE,IM,KIND=RK)
	MU = NEW_MU

END SUBROUTINE SETMU

!Setup for |GL/GT|.
subroutine setGLGT(GLGT_in)
    implicit none
    real(kind=rk),intent(in) :: GLGT_in
    GLGT = GLGT_in
end subroutine setGLGT

!Setup for thetaLT.
subroutine setTHETALT(THETALT_in)
    implicit none
    real(kind=rk),intent(in) :: THETALT_in
    THETALT = THETALT_in
end subroutine setTHETALT

!Setup for form factors in the 3-body decays of omega and phi.
subroutine setFormFact(FormFact_in)
    implicit none
    logical,intent(in) :: FormFact_in
    FormFact = FormFact_in
end subroutine setFormFact

!************************************************!
!Here start the routines for the external decays.
!************************************************!

!************************************************!
!                 Two body decays                !
! VM -> PS + PS                                  !
! VM -> PS + V                                   !
!************************************************!

!Calculates the 4-dimensional scalar product.
FUNCTION FOUR_PROD(P1,P2)
    IMPLICIT NONE
    REAL(KIND=RK),DIMENSION(4),INTENT(IN)::P1,P2
    REAL(KIND=RK) :: FOUR_PROD
    
    FOUR_PROD = P1(4)*P2(4)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)
END FUNCTION FOUR_PROD

!Calculates the transverse mass.
FUNCTION TRANSVERSE_MASS(PHAD)
    IMPLICIT NONE
    REAL(KIND=RK),DIMENSION(4),INTENT(IN) :: PHAD
    REAL(KIND=RK) :: TRANSVERSE_MASS
    REAL(KIND=RK) :: M2
    
    M2=FOUR_PROD(PHAD,PHAD)
    TRANSVERSE_MASS=SQRT(M2+PHAD(1)*PHAD(1)+PHAD(2)*PHAD(2))
END FUNCTION TRANSVERSE_MASS

!Boosts the 4-vector of the daughter from the rest frame of
!the mother to the frame where the mother was produced.
SUBROUTINE BOOST_MATRIX(PMOTHER,PDAUGHTER,PDAUGHTER_BOOSTED)
    IMPLICIT NONE
    REAL(KIND=RK),DIMENSION(4),INTENT(IN)::PMOTHER,PDAUGHTER
    REAL(KIND=RK),DIMENSION(4),INTENT(OUT)::PDAUGHTER_BOOSTED
    REAL(KIND=RK) :: MM,MMT,MPX,MPY,MPT2
    REAL(KIND=RK) :: T00,T01,T02,T10,T11,T12,T20,T21,T22
    REAL(KIND=RK) :: PD0,PD1,PD2,PD3
    REAL(KIND=RK) :: L00,L03,L30,L33
    
    !In 4-vectors the time component is the last component.
    MM = SQRT(FOUR_PROD(PMOTHER,PMOTHER))
    MMT = TRANSVERSE_MASS(PMOTHER)
    MPX = PMOTHER(1)
    MPY = PMOTHER(2)
    MPT2 = MPX*MPX+MPY*MPY
    
    !Transverse boost.
    T00 = MMT/MM
    T01 = MPX/MM
    T02 = MPY/MM
    T10 = MPX/MM
    T11 = 1._RK+(MMT/MM - 1._RK)*(MPX*MPX/MPT2)
    T12 = (MMT/MM-1._RK)*(MPX*MPY/MPT2) 
    T20 = MPY/MM
    T21 = T12
    T22 = 1._RK+(MMT/MM - 1._RK)*(MPY*MPY/MPT2)
    !New momentum of the daughter.
    PD0 = T00*PDAUGHTER(4)+T01*PDAUGHTER(1)+T02*PDAUGHTER(2)
    PD1 = T10*PDAUGHTER(4)+T11*PDAUGHTER(1)+T12*PDAUGHTER(2)
    PD2 = T20*PDAUGHTER(4)+T21*PDAUGHTER(1)+T22*PDAUGHTER(2)
    PD3 = PDAUGHTER(3)
    
    !Longitudinal boost.
    L00 = PMOTHER(4)/MMT
    L03 = PMOTHER(3)/MMT
    L30 = L03
    L33 = L00
    !Final momentum of the daughter.
    PDAUGHTER_BOOSTED(4) = L00*PD0+L03*PD3
    PDAUGHTER_BOOSTED(1) = PD1
    PDAUGHTER_BOOSTED(2) = PD2
    PDAUGHTER_BOOSTED(3) = L30*PD0+L33*PD3

END SUBROUTINE BOOST_MATRIX

!Subroutine to calculate the density matrix of a vector meson
!in the splitting q->VM+q' according to the rules of the string+3P0 model.
subroutine VMdensityMatrix(Sx,Sy,Sz,kpx,kpy,SDM_had)
    implicit none
    real(kind=rk),intent(in) :: Sx,Sy,Sz,kpx,kpy
    type(SpinDensityMatrix),intent(out) :: SDM_had
    !Matrices for the calculation of the density matrix of q and 3P0 propagator.
    complex(kind=rk),dimension(2,2) :: RHO, DELTA, DELTADAG
    !Matrices for the calculation of the VM spin density matrix
    COMPLEX(kind=rk), DIMENSION(2,2) :: XRHOX,XRHOY,XRHOZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: YRHOX,YRHOY,YRHOZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: ZRHOX,ZRHOY,ZRHOZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: MRHOXX,MRHOXY,MRHOXZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: MRHOYX,MRHOYY,MRHOYZ
    COMPLEX(kind=rk), DIMENSION(2,2) :: MRHOZX,MRHOZY,MRHOZZ
    COMPLEX(kind=rk)                 :: RHOXX,RHOXY,RHOXZ
    COMPLEX(kind=rk)                 :: RHOYX,RHOYY,RHOYZ
    COMPLEX(kind=rk)                 :: RHOZX,RHOZY,RHOZZ
    COMPLEX(kind=rk)                 :: TRACE
    real(kind=rk),dimension(3,3)     :: rhoHad

    !Spin density matrix of q.
    rho = 0.5 * (iden + Sx * sig1 &
                      + Sy * sig2 &
                      + Sz * sig3 )

    !Build up the 3P0 propagator.
    DELTA = mu*IDEN + matmul(SIG3,SIG1)*kpx + matmul(SIG3,SIG2)*kpy
    DELTADAG = conjg(mu)*IDEN - matmul(SIG3,SIG1)*kpx - matmul(SIG3,SIG2)*kpy

    !Calculate the density matrix of the VM.
    !Multiply by the VM emission vertex (coupling constants added later).
    XRHOX = matmul(SIG1, matmul(RHO,SIG1))
    XRHOY = matmul(SIG1, matmul(RHO,SIG2))
    XRHOZ = matmul(SIG1, RHO)

    YRHOY = matmul(SIG2, matmul(RHO,SIG2))
    YRHOX = matmul(SIG2, matmul(RHO,SIG1))
    YRHOZ = matmul(SIG2, RHO)

    ZRHOZ = RHO
    ZRHOX = matmul(RHO, SIG1)
    ZRHOY = matmul(RHO, SIG2)

    !Multiply by the 3P0 propagator and coupling constants.
    MRHOXX =  matmul(DELTA, matmul(YRHOY,DELTADAG))
    MRHOXY = -matmul(DELTA, matmul(YRHOX,DELTADAG))
    MRHOXZ = -matmul(DELTA, matmul(YRHOZ,DELTADAG))* (GLGT) * (SIN(THETALT)+IM*COS(THETALT))

    MRHOYY =  matmul(DELTA, matmul(XRHOX,DELTADAG))
    MRHOYX = -matmul(DELTA, matmul(XRHOY,DELTADAG))
    MRHOYZ =  matmul(DELTA, matmul(XRHOZ,DELTADAG)) * (GLGT) * (SIN(THETALT)+IM*COS(THETALT))

    MRHOZZ = matmul(DELTA, matmul(ZRHOZ,DELTADAG)) * (GLGT**2)
    MRHOZX = matmul(DELTA, matmul(ZRHOY,DELTADAG)) * (GLGT)*(-SIN(THETALT)+IM*COS(THETALT))
    MRHOZY = matmul(DELTA, matmul(ZRHOX,DELTADAG)) * (GLGT)*(SIN(THETALT)-IM*COS(THETALT))

    !VM spin density matrix elements in the rest frame.
    RHOXX = MRHOXX(1,1)+MRHOXX(2,2)
    RHOXY = MRHOXY(1,1)+MRHOXY(2,2)
    RHOXZ = MRHOXZ(1,1)+MRHOXZ(2,2)

    RHOYY = MRHOYY(1,1)+MRHOYY(2,2)
    RHOYX = MRHOYX(1,1)+MRHOYX(2,2)
    RHOYZ = MRHOYZ(1,1)+MRHOYZ(2,2)

    RHOZZ = MRHOZZ(1,1)+MRHOZZ(2,2)
    RHOZX = MRHOZX(1,1)+MRHOZX(2,2)
    RHOZY = MRHOZY(1,1)+MRHOZY(2,2)

    !Normalize the density matrix dividing by the trace.
    TRACE = RHOXX+RHOYY+RHOZZ

    RHOXX = RHOXX/TRACE
    RHOXY = RHOXY/TRACE
    RHOXZ = RHOXZ/TRACE
    RHOYX = RHOYX/TRACE
    RHOYY = RHOYY/TRACE
    RHOYZ = RHOYZ/TRACE
    RHOZX = RHOZX/TRACE
    RHOZY = RHOZY/TRACE
    RHOZZ = RHOZZ/TRACE

    !For the considered decays only the symmetric part of the density matrix is needed.
    !It coincides with the real part of the matrix.
    rhoHad(1,1) = real(RHOXX)
    rhoHad(1,2) = real(RHOXY)
    rhoHad(1,3) = real(RHOXZ)
    rhoHad(2,1) = real(RHOYX)
    rhoHad(2,2) = real(RHOYY)
    rhoHad(2,3) = real(RHOYZ)
    rhoHad(3,1) = real(RHOZX)
    rhoHad(3,2) = real(RHOZY)
    rhoHad(3,3) = real(RHOZZ)
    !Set up the density matrix of the hadron.
    call SDM_had%initElements( rhoHad(1,1), rhoHad(1,2), rhoHad(1,3),&
                               rhoHad(2,1), rhoHad(2,2), rhoHad(2,3),&
                               rhoHad(3,1), rhoHad(3,2), rhoHad(3,3)  )

end subroutine VMdensityMatrix

!Gives the energy and the absolute momentum of the daughters
!in a 2-body decay in the mother's rest frame.
subroutine TwoBodyDecay(MM,M1,M2,E1,E2,P)
    implicit none
    real(kind=rk),intent(in) :: MM,M1,M2
    real(kind=rk),intent(out):: E1,E2,P
    REAL(KIND=RK) :: MS,MD

    MS = M1 + M2
    MD = M1 - M2
    P  = sqrt((MM*MM-MD*MD) * (MM*MM-MS*MS))/(2.0*MM)
    E1 = sqrt(P*P + M1*M1)
    E2 = sqrt(P*P + M2*M2)
end subroutine TwoBodyDecay

!Subroutine for the generation of decay_theta, the polar angle
!of the relative momentum of the decay products in a 2-body decay.
subroutine Gen_Decay_Theta(rhozz,decay_theta)
    implicit none
    real(kind=rk),intent(in) :: rhozz
    real(kind=rk),intent(out):: decay_theta
    real(kind=rk) :: rnd, rnd2, asym, xx, w, wL, wC, wR
    
    if(rhozz>1.) print*,"mc3P0:: Gen_Decay_Theta in infinite loop. TYPE CTRL+C."
    asym = -(1._rk - 3._rk * rhozz)/(1._rk - rhozz)

    do
        call random_number(w)
        call random_number(rnd)
        call random_number(rnd2)
        !Case asym > 0.
        if( asym > 0. ) then
            if( w > 0.5) then
                xx = (-1._rk + sqrt(1._rk + asym*(2._rk+asym)*rnd))/asym
                rnd2 = rnd2 *(1._rk + asym*xx)
            else
                xx = (1._rk - sqrt(1._rk + asym*(2._rk+asym)*(1._rk-rnd)))/asym
                rnd2 = rnd2 *(1._rk - asym*xx)
            endif
        !Case asym < 0
        else
            wC = 0.25_rk*(12. + asym)/(6. + 2.*asym)
            wL = (1. - wC)/2
            wR = wL
            if( w < wC ) then
                xx = -1._rk/2 + rnd
                rnd2 = rnd2 * 1._rk
            else if( w > wC .and. w < (wC+wR) ) then
                xx = asym - 1._rk + sqrt(1._rk + asym*(2._rk + asym)*rnd)
                xx = 0.5_rk * xx/asym
                rnd2 = rnd2 * (2._rk*asym*xx + 1._rk - asym)
            else
                xx = 1._rk - asym - sqrt((1._rk + asym)**2 - asym*(2._rk + asym)*rnd)
                xx = 0.5_rk * xx/asym
                rnd2 = rnd2 * (-2._rk*asym*xx + 1._rk - asym)
            endif
        endif
            if( rnd2 < (1._rk + asym * xx**2)) then
                decay_theta = acos(xx)
                if((1._rk-rhozz-(1._rk-3._rk*rhozz)*xx*xx)<0) then
                    print*,"mc3P0:: Negative distribution in Gen_Decay_Theta!"
                    print*,"ABORTING."
                    stop
                endif
                exit
            endif
    enddo
end subroutine Gen_Decay_Theta

!Function to evaluate the probability density of the azimuthal angle
!phi in a two-body decay.
FUNCTION ANGULAR_PHI(XROW,YROW,ZROW,DECAY_THETA,DECAY_PHI)
    IMPLICIT NONE
    REAL(KIND=RK),DIMENSION(3),INTENT(IN) :: XROW,YROW,ZROW
    REAL(KIND=RK), INTENT(IN) :: DECAY_THETA,DECAY_PHI
    REAL(KIND=RK) :: ANGULAR_PHI
    REAL(KIND=RK) :: SIN2,COS2
    REAL(KIND=RK) :: W0,W1,W2,W3,W4

    SIN2=SIN(DECAY_THETA)*SIN(DECAY_THETA)
    COS2=COS(DECAY_THETA)*COS(DECAY_THETA)
    
    W0=YROW(2)*SIN2+ZROW(3)*COS2
    W1=(XROW(1)-YROW(2))*SIN2
    W2=XROW(3)*SIN(2._RK*DECAY_THETA)
    W3=YROW(3)*SIN(2._RK*DECAY_THETA)
    W4=XROW(2)*SIN2
    
    ANGULAR_PHI=W0+W1*(COS(DECAY_PHI)**2)+W2*COS(DECAY_PHI)+W3*SIN(DECAY_PHI)+W4*SIN(2*DECAY_PHI)
END FUNCTION ANGULAR_PHI

!Generates the azimuthal angle phi in a two-body decay.
SUBROUTINE GEN_DECAY_PHI(XROW,YROW,ZROW,DECAY_THETA,DECAY_PHI)
    IMPLICIT NONE
    REAL(KIND=RK),DIMENSION(3),INTENT(IN) :: XROW,YROW,ZROW
    REAL(KIND=RK), INTENT(IN) :: DECAY_THETA
    REAL(KIND=RK), INTENT(OUT) :: DECAY_PHI
    REAL(KIND=RK) ::RND,RND2,SIN2,COS2
    REAL(KIND=RK) :: W0,W1,W2,W3,W4,MAXX

    SIN2=SIN(DECAY_THETA)*SIN(DECAY_THETA)
    COS2=COS(DECAY_THETA)*COS(DECAY_THETA)
    
    W0=YROW(2)*SIN2+ZROW(3)*COS2
    W1=ABS((XROW(1)-YROW(2))*SIN2)
    W2=ABS(XROW(3)*SIN(2._RK*DECAY_THETA))
    W3=ABS(YROW(3)*SIN(2._RK*DECAY_THETA))
    W4=ABS(XROW(2)*SIN2)
    
    MAXX=W0+W1+W2+W3+W4
    
    DO 
        CALL RANDOM_NUMBER(RND)
        RND=RND*2._RK*PI
        CALL RANDOM_NUMBER(RND2)
        RND2=RND2*MAXX
        IF(RND2<ANGULAR_PHI(XROW,YROW,ZROW,DECAY_THETA,RND)) THEN
            DECAY_PHI=RND
            EXIT
        ENDIF
    ENDDO
END SUBROUTINE GEN_DECAY_PHI

!Generates a complete two-body decay of the type VM->PS+PS
!in the VM rest frame.
subroutine GenVtoPsPsDecay3D_new(Mres,M1,M2,rhoIn,p1Vec4,relVec)
    implicit none
    real(kind=rk),intent(in)                :: Mres,M1,M2
    type(SpinDensityMatrix) :: rhoIn
    real(kind=rk),intent(out),dimension(4)  :: p1Vec4
    real(kind=rk),intent(out),dimension(3)  :: relVec
    real(kind=rk),dimension(3)              :: xRow,yRow,zRow
    real(kind=rk)                           :: decayTheta,decayPhi
    real(kind=rk)                           :: E1,E2,pMod

    ! Save thw rows of the spin density matrix.
    !x-row.
    xRow(1)=REAL(rhoIn%xx)
    xRow(2)=REAL(rhoIn%xy)
    xRow(3)=REAL(rhoIn%xz)
    !y-row.
    yRow(1)=REAL(rhoIn%yx)
    yRow(2)=REAL(rhoIn%yy)
    yRow(3)=REAL(rhoIn%yz)
    !z-row.
    zRow(1)=REAL(rhoIn%zx)
    zRow(2)=REAL(rhoIn%zy)
    zRow(3)=REAL(rhoIn%zz)

    !Generate energies and momenta of the daughters in the mother's rest frame.
    call TwoBodyDecay(Mres,M1,M2,E1,E2,pMod)

    !Generate the decay polar angle whose distribution depends only on rhozz.
    call GEN_DECAY_THETA(ZROW(3),decayTheta)

    !Generate the decay azimuthal angle whose distribution depends
    !on all density matrix elements.
    call GEN_DECAY_PHI(XROW,YROW,ZROW,decayTheta,decayPhi)

    !Relative momentum of the daughters in the mother's rest frame.
    relVec(1) = COS(decayPhi) * SIN(decayTheta)
    relVec(2) = SIN(decayPhi) * SIN(decayTheta)
    relVec(3) = COS(decayTheta)

    !Four momentum of daughter 1 in the mother's rest frame.
    !(that of daughter 2 can be calculated by energy-momentum conservation).
    p1Vec4(1) = pMod*relVec(1)
    p1Vec4(2) = pMod*relVec(2)
    p1Vec4(3) = pMod*relVec(3)
    p1Vec4(4) = E1
end subroutine GenVtoPsPsDecay3D_new

!Subroutine for the generation of the VM -> PS + PS decay.
subroutine VMtoPSPSdecay(Px,Py,Pz,E,&
                         M,M1,M2,&
                         Sx,Sy,Sz,&
                         kpx,kpy,&
                         P1x,P1y,P1z,E1,&
                         P2x,P2y,P2z,E2,&
                         Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,&
                         Dzx,Dzy,Dzz)
    implicit none
    real(kind=rk),intent(in) :: Px,Py,Pz,E,M,M1,M2
    real(kind=rk),intent(in) :: Sx,Sy,Sz
    real(kind=rk),intent(in) :: kpx,kpy
    real(kind=rk),intent(out):: P1x,P1y,P1z,E1
    real(kind=rk),intent(out):: P2x,P2y,P2z,E2
    real(kind=rk),intent(out):: Dxx,Dxy,Dxz
    real(kind=rk),intent(out):: Dyx,Dyy,Dyz
    real(kind=rk),intent(out):: Dzx,Dzy,Dzz
    !Variables for internal use.
    type(SpinDensityMatrix)  :: rhoHad
    real(kind=rk),dimension(3,3) :: rhoCheck
    real(kind=rk),dimension(4)   :: p1, p2, pRes
    real(kind=rk),dimension(3)   :: relVec

    !Calculate the VM density matrix.
    call VMdensityMatrix(Sx,Sy,Sz,kpx,kpy,rhoHad)

    !Generate the momenta of the daughters in the mother's rest frame.
    call GenVtoPsPsDecay3D_new(M,M1,M2,rhoHad,p1,relVec)

    !Boost the daughters to the string frame.
    pRes(1) = Px
    pRes(2) = Py
    pRes(3) = Pz
    pRes(4) = E
    CALL BOOST_MATRIX(pRes, p1, p1)
    p2 = pRes - p1

    !4-momentum of daughter 1.
    P1x = p1(1)
    P1y = p1(2)
    P1z = p1(3)
    E1  = p1(4)
    !4-momentum of daughter 2.
    P2x = p2(1)
    P2y = p2(2)
    P2z = p2(3)
    E2  = p2(4)

    !Check that energy calculation is ok.
    if( (p1(4)-pRes(4)) .gt. 0.0001 ) PRINT*,"mc3P0:: ERROR in VM->PS+PS: Edaughter1 > Emother."
    if( (p2(4)-pRes(4)) .gt. 0.0001 ) PRINT*,"mc3P0:: ERROR in VM->PS+PS: Edaughter2 > Emother."

    !Acceptance density matrix.
    Dxx = relVec(1) * relVec(1)
    Dxy = relVec(1) * relVec(2)
    Dxz = relVec(1) * relVec(3)
    Dyx = relVec(2) * relVec(1)
    Dyy = relVec(2) * relVec(2)
    Dyz = relVec(2) * relVec(3)
    Dzx = relVec(3) * relVec(1)
    Dzy = relVec(3) * relVec(2)
    Dzz = relVec(3) * relVec(3)
end subroutine VMtoPSPSdecay

!Subroutine for the generation of the VM -> PS + V decay.
subroutine VMtoPSVdecay(Px,Py,Pz,E,&
                         M,M1,M2,&
                         Sx,Sy,Sz,&
                         kpx,kpy,&
                         P1x,P1y,P1z,E1,&
                         P2x,P2y,P2z,E2,&
                         Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,&
                         Dzx,Dzy,Dzz)
    implicit none
    real(kind=rk),intent(in) :: Px,Py,Pz,E,M,M1,M2
    real(kind=rk),intent(in) :: Sx,Sy,Sz
    real(kind=rk),intent(in) :: kpx,kpy
    real(kind=rk),intent(out):: P1x,P1y,P1z,E1
    real(kind=rk),intent(out):: P2x,P2y,P2z,E2
    real(kind=rk),intent(out):: Dxx,Dxy,Dxz
    real(kind=rk),intent(out):: Dyx,Dyy,Dyz
    real(kind=rk),intent(out):: Dzx,Dzy,Dzz
    !Variables for internal use.
    type(SpinDensityMatrix)  :: rhoHad
    real(kind=rk),dimension(3,3) :: rhoCheck
    real(kind=rk),dimension(4)   :: p1, p2, pRes
    real(kind=rk),dimension(3)   :: relVec

    !Calculate the VM density matrix.
    call VMdensityMatrix(Sx,Sy,Sz,kpx,kpy,rhoHad)

    !Generate the momenta of the daughters in the mother's rest frame.
    call GenVtoPsPsDecay3D_new(M,M1,M2,rhoHad,p1,relVec)

    !Boost the daughters to the string frame.
    pRes(1) = Px
    pRes(2) = Py
    pRes(3) = Pz
    pRes(4) = E
    CALL BOOST_MATRIX(pRes, p1, p1)
    p2 = pRes - p1

    !4-momentum of daughter 1.
    P1x = p1(1)
    P1y = p1(2)
    P1z = p1(3)
    E1  = p1(4)

    !4-momentum of daughter 2.
    P2x = p2(1)
    P2y = p2(2)
    P2z = p2(3)
    E2  = p2(4)

    !Check that energy calculation is ok.
    if( (p1(4)-pRes(4)) .gt. 0.0001 ) PRINT*,"mc3P0:: ERROR in VM->V+PS: Edaughter1 > Emother."
    if( (p2(4)-pRes(4)) .gt. 0.0001 ) PRINT*,"mc3P0:: ERROR in VM->V+PS: Edaughter2 > Emother."
    
    !Acceptance density matrix.
    Dxx =  1.0_rk - relVec(1 ) * relVec(1)
    Dxy = -relVec(1) * relVec(2)
    Dxz = -relVec(1) * relVec(3)
    Dyx = -relVec(2) * relVec(1)
    Dyy =  1.0_rk - relVec(2) * relVec(2)
    Dyz = -relVec(2) * relVec(3)
    Dzx = -relVec(3) * relVec(1)
    Dzy = -relVec(3) * relVec(2)
    Dzz =  1.0_rk - relVec(3) * relVec(3)
end subroutine VMtoPSVdecay

!************************************************!
!                 Three body decays              !
! omega/phi -> PS + PS + PS                      !
!************************************************!

!Generate energies and momenta of the three decay products.
subroutine ThreeEnergyMomentum(Mres,E1,E2,E3,p1,p2,p3,m12,m13,m23)
    implicit none
    real(kind=rk),intent(in)  :: Mres
    real(kind=rk),intent(out) :: E1,E2,E3,p1,p2,p3,m12,m13,m23
    real(kind=rk)              :: m1,m2,m3
    real(kind=rk)              :: x1,x2,x3
    real(kind=rk)              :: x1Min,x1Max,x2Min,x2Max,x3Min,x3Max
    real(kind=rk)              :: mAng, y
    real(kind=rk)              :: p1Trial,p2Trial,p3Trial,cosTheta23
    real(kind=rk)              :: E1max,E2max,E3max

    m1 = MPI    !pi+
    m2 = MPI    !pi-
    m3 = MPIZ   !pi0
    
    E1max = (Mres**2 + Mpi**2 - (Mpi+Mpiz)**2)/(2.0*Mres)
    E2max = (Mres**2 + Mpi**2 - (Mpi+Mpiz)**2)/(2.0*Mres)
    E3max = (Mres**2 + Mpiz**2 - (2*Mpi)**2)/(2.0*Mres)
    
    x1Min = m1/Mres
    x1Max = E1max/Mres
    
    x2Min = m2/Mres
    x2Max = E2max/Mres
    
    x3Min = m3/Mres
    x3Max = E3max/Mres
    
    do
        !Energy conservation by generating the triangle
        !with sides x1, x2 and x3.
        call random_number(x1)
        x1 = x1Min + x1 * (x1Max - x1Min)
        call random_number(x2)
        x2 = x2Min + x2 * (x2Max - x2Min)
        if(x2 < (0.5_rk - x1)) then
            x1 = x1Max - x1
            x2 = x2Max - x2
        endif
        x3 = 1.0_rk - x1 - x2

        !Impose momentum conservation.
        if ( x1*x1*Mres*Mres - m1*m1 < 0.0 ) cycle
        if ( x2*x2*Mres*Mres - m2*m2 < 0.0 ) cycle
        if ( x3*x3*Mres*Mres - m3*m3 < 0.0 ) cycle
        p1Trial = sqrt(x1*x1*Mres*Mres - m1*m1)
        p2Trial = sqrt(x2*x2*Mres*Mres - m2*m2)
        p3Trial = sqrt(x3*x3*Mres*Mres - m3*m3)
        
        cosTheta23 = p1Trial*p1Trial - p2Trial*p2Trial - p3Trial*p3Trial
        cosTheta23 = cosTheta23 / (2.0 * p2Trial * p3Trial)
        
        if( abs(cosTheta23) <= 1.0 .and. (x3 > x3Min) .and. (x3 < x3Max) ) then
            E1 = x1*Mres
            E2 = x2*Mres
            E3 = x3*Mres
            
            p1 = sqrt(E1*E1 - m1*m1)
            p2 = sqrt(E2*E2 - m2*m2)
            p3 = sqrt(E3*E3 - m3*m3)
            
            m12 = sqrt(Mres*Mres + m3*m3 - 2.0*Mres*E3)
            m13 = sqrt(Mres*Mres + m2*m2 - 2.0*Mres*E2)
            m23 = sqrt(Mres*Mres + m1*m1 - 2.0*Mres*E1)
            exit
        endif
    enddo
end subroutine ThreeEnergyMomentum

! The matrix element squared for a three-body decay.
function ThreeMatrixElement(p1,p2,p3)
    implicit none
    real(kind=rk)                ::ThreeMatrixElement
    real(kind=rk),intent(in)    :: p1,p2,p3
    real(kind=rk)                :: cosTheta13

    cosTheta13 = (p2*p2 - p3*p3 - p1*p1) / (2.0 * p3 * p1)
    ThreeMatrixElement = p3 * p3 * p1 * p1 * (1.0_rk - cosTheta13 * cosTheta13)
end function ThreeMatrixElement

!The Breit-Wigner function for a rho meson.
!Needed only in three body decays.
function BreitWignerFunc(Mres,KFres,Gammares)
    implicit none
    real(kind=rk)                 :: BreitWignerFunc
    real(kind=rk),intent(in)     :: Mres,Gammares
    integer,intent(in)            :: KFres

    if(KFres == KFrhop .OR. KFres == KFrhom) then
        BreitWignerFunc = 1.0_rk / ((Mres*Mres - Mrho*Mrho)**2 + (Mrho*Gammares)**2)
    else if(KFres == KFrhoz) then
        BreitWignerFunc = 1.0_rk / ((Mres*Mres - Mrhoz*Mrhoz)**2 + (Mrhoz*Gammares)**2)
    else
        print*,"Error in mc3P0::BreittWignerFunc."
        print*,"Intermediate resonance not implemented."
        print*,"Stop."
        stop
    endif
end function BreitWignerFunc

!The form factor |BW1+BW2+BW3|^2 where 1,2,3 are the intermediate
!rho mesons in the three-body decays of omega and phi.
function ThreeBreitWigner(m12,m13,m23)
    implicit none
    real(kind=rk)                :: ThreeBreitWigner
    real(kind=rk),intent(in)    :: m12,m13,m23
    real(kind=rk)                :: m0,mp,mm
    real(kind=rk)                :: qz,qpm,mSum,mDiff
    real(kind=rk)                :: gRhozPiPi,gRhopPiPi,gRhomPiPi
    real(kind=rk)                :: szz,spp,smm,szp,szm,spm
    real(kind=rk)                :: term1,term2,term3
    
    m0 = m12
    mp = m13
    mm = m23
    
    qz = 0.5 * sqrt(Mrhoz*Mrhoz - 4.0*Mpi*Mpi)
    mSum = Mpiz + Mpi
    mDiff = Mpi - Mpiz
    qpm  = 0.5 * sqrt((Mrho*Mrho - mSum*mSum) * (Mrho*Mrho - mDiff*mDiff)) / Mrho
    
    gRhozPiPi = sqrt( 6.0 * pi * Mrhoz * Mrhoz * GammaRhozPiPi / (qz * qz * qz) )
    gRhopPiPi  = sqrt( 6.0 * pi * Mrho * Mrho * GammaRhopPiPi / (qpm * qpm * qpm) )
    gRhomPiPi  = sqrt( 6.0 * pi * Mrho * Mrho * GammaRhomPiPi / (qpm * qpm * qpm) )
    
    szz = BreitWignerFunc(m0,KFrhoz,GammaRhozPiPi)
    spp = BreitWignerFunc(mp,KFrhop,GammaRhopPiPi)
    smm = BreitWignerFunc(mm,KFrhom,GammaRhomPiPi)
    szp = szz * spp * ((m0**2 - Mrhoz**2) * (mp**2 - Mrho**2) - Mrhoz*GammaRhozPiPi*Mrho*GammaRhopPiPi)
    szm = szz * smm * ((m0**2 - Mrhoz**2) * (mm**2 - Mrho**2) - Mrhoz*GammaRhozPiPi*Mrho*GammaRhomPiPi)
    spm = spp * smm * ((mp**2 - Mrho**2) * (mm**2 - Mrho**2) - Mrho*GammaRhopPiPi*Mrho*GammaRhomPiPi)
    
    term1 = gRhozPiPi * gRhozPiPi * szz + gRhopPiPi * gRhopPiPi * spp + gRhomPiPi * gRhomPiPi * smm
    term2 = 2.0 * gRhozPiPi * (gRhopPiPi * szp + gRhomPiPi * szm)
    term3 = 2.0 * gRhopPiPi * gRhomPiPI * spm
    
    ThreeBreitWigner =  term1 + term2 + term3
end function ThreeBreitWigner

! Total weight (matrix element x form factor) for a three-body decay.
function ThreeBodyWeight(Mres,p1,p2,p3,m12,m13,m23)
    implicit none
    real(kind=rk)                :: ThreeBodyWeight
    real(kind=rk),intent(in)    :: Mres,p1,p2,p3,m12,m13,m23
    real(kind=rk)                :: Gamma0,p1max2,p3max2,xx,yy,zz
    real(kind=rk)                :: m0,mp,mm
    real(kind=rk)                :: qz,qpm,mSum,mDiff
    real(kind=rk)                :: gRhozPiPi,gRhopPiPi,gRhomPiPi
    real(kind=rk)                :: MatEMax,maxpp,maxmm,maxzz,FormFactMax,delta,p2Equi
    
    !Calculate the maximum value of the matrix element.
    xx = Mres**2 - (2*Mpi+Mpiz)**2
    yy = Mres**2 - Mpiz**2
    zz = Mres**2 - (Mpiz-2*Mpi)**2
    
    p1max2 = xx * yy / (4*Mres*Mres)
    p3max2 = xx * zz / (4*Mres*Mres)
    MatEMax = p1max2 * p3max2
    
    !Calculate the maximum value of the form factor.
    if(FormFact) then
        m0 = m12
        mp = m13
        mm = m23
    
        qz = 0.5 * sqrt(Mrhoz*Mrhoz - 4.0*Mpi*Mpi)
        mSum = Mpiz + Mpi
        mDiff = Mpi - Mpiz
        qpm  = 0.5 * sqrt((Mrho*Mrho - mSum*mSum) * (Mrho*Mrho - mDiff*mDiff)) / Mrho
    
        gRhozPiPi = sqrt( 6.0 * pi * Mrhoz * Mrhoz * GammaRhozPiPi / (qz * qz * qz) )
        gRhopPiPi  = sqrt( 6.0 * pi * Mrho * Mrho * GammaRhopPiPi / (qpm * qpm * qpm) )
        gRhomPiPi  = sqrt( 6.0 * pi * Mrho * Mrho * GammaRhomPiPi / (qpm * qpm * qpm) )
    
        maxpp = gRhopPiPi/(Mrho*GammaRhopPiPi)
        maxmm = gRhomPiPi/(Mrho*GammaRhomPiPi)
        maxzz = gRhozPiPi/(Mrhoz*GammaRhozPiPi)
    
        FormFactMax = (maxpp + maxmm + maxzz)**2
    endif    
    !Calculate the total three-body weight = ME x FF / max(ME) max(FF).
    !If FormFact=false, the form factor is not included in the weight.
    if(FormFact) then
        ThreeBodyWeight = ThreeMatrixElement(p1,p2,p3) * ThreeBreitWigner(m12,m13,m23)
        ThreeBodyWeight = ThreeBodyWeight / ( MatEMax * FormFactMax )
    else
        ThreeBodyWeight = ThreeMatrixElement(p1,p2,p3)
        ThreeBodyWeight = ThreeBodyWeight / MatEMax
    endif
    
    if(ThreeBodyWeight>1) then
        print*,"mc3P0::ThreeBodyWeight: Error in calculation of the weight."
        if(FormFact) then
            print*,"ThreeBreitWigner: ",ThreeBreitWigner(m12,m13,m23)
            print*,"FormFactMax: ",FormFactMax
        endif
        print*,"ThreeBodyWeight: ",ThreeBodyWeight
    endif
    
    !Error message.
    if( p1>sqrt( p1max2 ) ) print*,"mc3P0::ThreeBodyWeight: p1>p1max."
    if( p3>sqrt( p3max2 ) ) print*,"mc3P0::ThreeBodyWeight: p3>p3max."
end function ThreeBodyWeight


!Generates invariant masses for a VM -> pi pi pi decay.
subroutine GenThreeBodyDecay_new(Mres,&
                                 E1,E2,E3,&
                                 p1,p2,p3,&
                                 m12,m13,m23)
    implicit none
    real(kind=rk),intent(in)  :: Mres
    real(kind=rk),intent(out) :: E1,E2,E3,m12,m13,m23,p1,p2,p3
    real(kind=rk)             :: xm12,xm13,xm23,xp1,xp2,xp3
    real(kind=rk)             :: rnd1,rnd2,rnd3,rnd
    real(kind=rk)             :: max,weight

    do
      !Generate energies and momenta according to a Dalitz decay.
      call ThreeEnergyMomentum(Mres,rnd1,rnd2,rnd3,xp1,xp2,xp3,xm12,xm13,xm23)

      !Calculate weight according to form factors for phi/omega -> 3pi.
      weight = ThreeBodyWeight(Mres,xp1,xp2,xp3,xm12,xm13,xm23)
      call random_number(rnd)
      max = 1.0
      rnd = rnd * max
      if (weight > max) then
        print*,"mc3P0::GenThreeBodyDecay: weight larger than maximum."
      endif
      if(rnd < weight) then
        E1 = rnd1
        E2 = rnd2
        E3 = rnd3
        p1 = xp1
        p2 = xp2
        p3 = xp3
        m12 = xm12
        m13 = xm13
        m23 = xm23
        exit
      endif
    enddo

end subroutine GenThreeBodyDecay_new

!Generates the triangle of three-momenta for a Dalitz decay
!in the resonance rest frame.
subroutine GenDalitzTriangle(p1,p2,p3,p1Vec,p2Vec,p3Vec)
implicit none
real(kind=rk),intent(in)                :: p1,p2,p3
real(kind=rk),intent(out),dimension(3)    :: p1Vec,p2Vec,p3Vec
real(kind=rk)                            :: cosTheta12,cosTheta13
real(kind=rk)                            :: phi1,phi2,phi3

    cosTheta13 = (p2*p2 - p1*p1 - p3*p3) / (2.0 * p1 * p3)
    cosTheta12 = (p3*p3 - p1*p1 - p2*p2) / (2.0 * p1 * p2)
    
    !Dalitz triangle when p1 is along x.
    p1Vec(1) = p1
    p1Vec(2) = 0.0
    p1Vec(3) = 0.0
    
    p2Vec(1) = p2 * cosTheta12
    p2Vec(2) = -p2 * sqrt(1.0 - cosTheta12 * cosTheta12)
    p2Vec(3) = 0.0
    
    p3Vec(1) = p3 * cosTheta13
    p3Vec(2) = p3 * sqrt(1.0 - cosTheta13 * cosTheta13)
    p3Vec(3) = 0.0
    
    !The azimuthal angle of p1 is random.
    call random_number(phi1)
    
    call RotzAxis(phi1,p1Vec,p1Vec)
    call RotzAxis(phi1,p2Vec,p2Vec)
    call RotzAxis(phi1,p3Vec,p3Vec)
end subroutine GenDalitzTriangle


!Rotation of an angle phi about z axis (for three-vectors).
subroutine RotzAxis(phi,VecIn,VecOut)
    implicit none
    real(kind=rk),intent(in)                :: phi
    real(kind=rk),intent(in),dimension(3)    :: VecIn
    real(kind=rk),intent(out),dimension(3)    :: VecOut
    real(kind=rk),dimension(3)                :: tempIn

    !Counterclockwise rotation of x towards y.
    tempIn = VecIn
    VecOut(1) = cos(phi) * tempIn(1) - sin(phi) * tempIn(2)
    VecOut(2) = sin(phi) * tempIn(1) + cos(phi) * tempIn(2)
    VecOut(3) = tempIn(3)
end subroutine RotzAxis


!Rotation of an angle theta about y axis (for three-vectors).
subroutine RotyAxis(theta,VecIn,VecOut)
    implicit none
    real(kind=rk),intent(in)                :: theta
    real(kind=rk),intent(in),dimension(3)    :: VecIn
    real(kind=rk),intent(out),dimension(3)    :: VecOut
    real(kind=rk),dimension(3)                :: tempIn

    !Counterclockwise rotation of z towards x.
    tempIn = VecIn
    VecOut(1) = cos(theta) * tempIn(1) + sin(theta) * tempIn(3)
    VecOut(2) = tempIn(2)
    VecOut(3) = -sin(theta) * tempIn(1) + cos(theta) * tempIn(3)
end subroutine RotyAxis


!Rotation of angles phi and theta (for three-vectors).
subroutine Rot3D(phi,theta,VecIn,VecOut)
    implicit none
    real(kind=rk),intent(in),dimension(3)    :: VecIn
    real(kind=rk),intent(in)                :: phi,theta
    real(kind=rk),intent(out),dimension(3)    :: VecOut

    !Rotate first about y of an angle theta.
    call RotyAxis(theta,VecIn,VecOut)

    !Then rotate about z of an angle phi.
    call RotzAxis(phi,VecIn,VecOut)
end subroutine Rot3D


!Generates the normal to a Dalitz triangle
subroutine GenNormalDalitzTriangle_new(rhoIn,NormalVec,phi,theta)
    implicit none
    type(SpinDensityMatrix)                 :: rhoIn
    real(kind=rk),intent(out),dimension(3)    :: NormalVec
    real(kind=rk),intent(out)                :: phi,theta
    real(kind=rk),dimension(3)                :: xRow,yRow,zRow

    !Save the rows of the spin density matrix.
    !x-row.
    xRow(1)=REAL(rhoIn%xx)
    xRow(2)=REAL(rhoIn%xy)
    xRow(3)=REAL(rhoIn%xz)
    !y-row.
    yRow(1)=REAL(rhoIn%yx)
    yRow(2)=REAL(rhoIn%yy)
    yRow(3)=REAL(rhoIn%yz)
    !z-row.
    zRow(1)=REAL(rhoIn%zx)
    zRow(2)=REAL(rhoIn%zy)
    zRow(3)=REAL(rhoIn%zz)

    !Generate the polar angle of the normal vector.
    call Gen_Decay_Theta( zRow(3),theta )
    
    !Generate the azimuthal angle of the normal vector.
    call Gen_Decay_Phi( xRow,yRow,zRow,theta,phi )

    !Construct the normal vector.
    NormalVec(1) = cos(phi) * sin(theta)
    NormalVec(2) = sin(phi) * sin(theta)
    NormalVec(3) = cos(theta)
end subroutine GenNormalDalitzTriangle_new


!Generates a three-body decay in the resonance rest frame.
subroutine GenThreeBodyDecay3D_new(Mres,rhoIn,&
                                   p1Vec4,p2Vec4,p3Vec4,&
                                   m12,m13,m23,&
                                   NormalVec)
implicit none
real(kind=rk),intent(in)                :: Mres
type(SpinDensityMatrix)                 :: rhoIn
real(kind=rk),intent(out),dimension(4)  :: p1Vec4,p2Vec4,p3Vec4
real(kind=rk),intent(out)               :: m12,m13,m23
real(kind=rk),intent(out),dimension(3)  :: NormalVec
real(kind=rk)                           :: E1,E2,E3,p1,p2,p3
real(kind=rk)                           :: phi,theta
real(kind=rk),dimension(3)              :: p1VecTemp,p2VecTemp,p3VecTemp
    
    !Generate energies, momenta and invariant masses according to phase space.
    call ThreeEnergyMomentum(Mres,E1,E2,E3,p1,p2,p3,m12,m13,m23)
    !Weight according to matrix element and form factor for omega/phi -> 3pi decay.
    call GenThreeBodyDecay_new(Mres,E1,E2,E3,p1,p2,p3,m12,m13,m23)
    !Generate the normal vector to the three-Body plane.
    call GenNormalDalitzTriangle_new(rhoIn,NormalVec,phi,theta)
    ! In the system where NormalVec is along z, generate the Dalitz triangle.
    call GenDalitzTriangle(p1,p2,p3,p1VecTemp,p2VecTemp,p3VecTemp)
    ! Rotate the decay momenta so that the normal to the plane is along NormalVec.
    call Rot3D(phi,theta,p1VecTemp,p1VecTemp)
    call Rot3D(phi,theta,p2VecTemp,p2VecTemp)
    call Rot3D(phi,theta,p3VecTemp,p3VecTemp)
    ! Final momenta in the resonance rest frame.
    p1Vec4(1) = p1VecTemp(1)
    p1Vec4(2) = p1VecTemp(2)
    p1Vec4(3) = p1VecTemp(3)
    p1Vec4(4) = E1
    
    p2Vec4(1) = p2VecTemp(1)
    p2Vec4(2) = p2VecTemp(2)
    p2Vec4(3) = p2VecTemp(3)
    p2Vec4(4) = E2
    
    p3Vec4(1) = p3VecTemp(1)
    p3Vec4(2) = p3VecTemp(2)
    p3Vec4(3) = p3VecTemp(3)
    p3Vec4(4) = E3

end subroutine GenThreeBodyDecay3D_new


!Subroutine for the generation of the VM -> PS + PS + PS decay
subroutine VMtoPSPSPSdecay(Px,Py,Pz,E,&
                         M,M1,M2,&
                         Sx,Sy,Sz,&
                         kpx,kpy,&
                         P1x,P1y,P1z,E1,&
                         P2x,P2y,P2z,E2,&
                         P3x,P3y,P3z,E3,&
                         Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,&
                         Dzx,Dzy,Dzz)
    implicit none
    real(kind=rk),intent(in) :: Px,Py,Pz,E,M,M1,M2
    real(kind=rk),intent(in) :: Sx,Sy,Sz
    real(kind=rk),intent(in) :: kpx,kpy
    real(kind=rk),intent(out):: P1x,P1y,P1z,E1
    real(kind=rk),intent(out):: P2x,P2y,P2z,E2
    real(kind=rk),intent(out):: P3x,P3y,P3z,E3
    real(kind=rk),intent(out):: Dxx,Dxy,Dxz
    real(kind=rk),intent(out):: Dyx,Dyy,Dyz
    real(kind=rk),intent(out):: Dzx,Dzy,Dzz
    !Variables for internal use.
    type(SpinDensityMatrix)  :: rhoHad
    real(kind=rk),dimension(3,3) :: rhoCheck
    real(kind=rk),dimension(4)   :: p1, p2, p3, pRes
    real(kind=rk),dimension(3)   :: NormalVec
    real(kind=rk) :: m12,m13,m23

    !Calculate the VM density matrix.
    call VMdensityMatrix(Sx,Sy,Sz,kpx,kpy,rhoHad)

    !Generate the momenta of the daughters in the mother's rest frame.
    call GenThreeBodyDecay3D_new(M,rhoHad,p1,p2,p3,m12,m13,m23,NormalVec)

    !Boost the daughters to the string frame.
    pRes(1) = Px
    pRes(2) = Py
    pRes(3) = Pz
    pRes(4) = E
    CALL BOOST_MATRIX(pRes,p1,p1)
    CALL BOOST_MATRIX(pRes,p2,p2)
    p3 = pRes - p1 - p2

    !4-momentum of daughter 1.
    P1x = p1(1)
    P1y = p1(2)
    P1z = p1(3)
    E1  = p1(4)
    !4-momentum of daughter 2.
    P2x = p2(1)
    P2y = p2(2)
    P2z = p2(3)
    E2  = p2(4)
    !4-momentum of daughter 3.
    P3x = p3(1)
    P3y = p3(2)
    P3z = p3(3)
    E3  = p3(4)

    !Check that energy calculation is ok.
    if( (p1(4)-pRes(4)) .gt. 0.0001 ) print*,"mc3P0::VMtoPsPsPsDecay: ERROR. Edaughter1 > Emother."
    if( (p2(4)-pRes(4)) .gt. 0.0001 ) print*,"mc3P0::VMtoPsPsPsDecay: ERROR. Edaughter2 > Emother."
    if( (p3(4)-pRes(4)) .gt. 0.0001 ) print*,"mc3P0::VMtoPsPsPsDecay: ERROR. Edaughter3 > Emother."

    !Acceptance density matrix.
    Dxx = NormalVec(1) * NormalVec(1)
    Dxy = NormalVec(1) * NormalVec(2)
    Dxz = NormalVec(1) * NormalVec(3)
    Dyx = NormalVec(2) * NormalVec(1)
    Dyy = NormalVec(2) * NormalVec(2)
    Dyz = NormalVec(2) * NormalVec(3)
    Dzx = NormalVec(3) * NormalVec(1)
    Dzy = NormalVec(3) * NormalVec(2)
    Dzz = NormalVec(3) * NormalVec(3)
end subroutine VMtoPSPSPSdecay


!End of module with routines for the implementation of the
!string+3P0 model.
END MODULE ROUTINES
