!
!
! TODO: 
!  - check sqrt(psi_h*psi_hbar) in the formulas, in many cases is cdabs(psi_h)
!  - check asymmetry factor in ML and PP (not appearing?)
!  - extinction depth is in intensity whereas pendellosung depth is in amplitude
!  - check the sign of c in ML Bragg
!
!



!----
!----
!---- MODULE:  crystal3
!----
!---- Main module for crystal diffraction
!----
!----
!----
!----
!----
!----      The sequence of calling is the following
!----
!----
!----      crystal_loadcrystaldata:  to load preprocessor file (now written by XOP)
!----      crystal_set:              to set (align) the crystal at a given energy
!----      crystal_info:             
!----      crysta_{prefect,mosaic,bent_pp,bent_ml} to compute reflectivity
!----
!----
!----      the interface is the single "crystal" routine that should be called with
!----      "kwhat" flag:
!----              <0 (calls crystal_loadcrystaldata)
!----               0 (calls crystal_set)
!----               1 (calls crystal_info)
!----               2 gets reflectivity
!----
!----      

Module crystal3
!---- Use Modules ----!

use shadow_globaldefinitions
use shadow_math
use elasticity

!---- Variables ----!
implicit none

!public

!todo move to a common place consistenly with shadow
real(kind=8) :: r0 = 2.817940289458D-13 ! classical electron radius in cm

!todo remove
!xianbos b_factor in c calculation (12.96 for ml model)(-0.0077 for pp model)
real(kind=8) :: B_factor

!
! type (structure) to store crystallographic data from file 
! (data usually created by a preprocessor)
!
integer(kind=ski),public,parameter        :: NMAXENER=1000, NMAXATOMS=100
type, public, bind(C) :: crystalData
    !C   RN: the constant (e^2/mc^2)/V or the ration between the classical e- 
    !C        radius and the volume of the unit cell [cm^(-2)]
    real(kind=skr)                               :: rn,d_spacing
    integer(kind=ski)                            :: itype
    real(kind=skr),dimension(NMAXATOMS)          :: atnum
    real(kind=skr),dimension(NMAXATOMS)          :: temper
    integer(kind=ski)                            :: nbatom
    complex(kind=skx),dimension(NMAXATOMS)       :: G,G_BAR
    real(kind=skr),dimension(NMAXATOMS)          :: G_0
    integer(kind=ski)                            :: nf0coeff
    real(kind=skr),dimension(NMAXATOMS,11)       :: f0coeff
    integer(kind=ski)                            :: npoint
    real(kind=skr),dimension(NMAXENER)           :: energy
    real(kind=skr),dimension(NMAXATOMS,NMAXENER) :: fp,fpp,fcompton
    real(kind=skr),dimension(NMAXATOMS)          :: fract
    character(len=sklen)                         :: filename
    real(kind=skr)                               :: version
end type crystalData


!
! this is the type that contains all the crystal information needed
! by crystal calculations. 
!
type, public, bind(C) :: oeSetup
    ! these are user inputs
    real(kind=skr)        :: thickness=0.1D0  ! crystal thickness 
    integer(kind=ski)     :: f_bragg_a=0      ! asymmetric?
    real(kind=skr)        :: a_bragg=0.0D0    ! asymmetry angle
    integer(kind=ski)     :: f_mosaic=0       ! [0] perfect, [1]mosaic, 
                                              ! [2] bent multilamellar
                                              ! [3] bent Penning-Polder
    real(kind=skr)        :: spread_mos=0.0D0 ! crystal mosaicity
    integer(kind=ski)     :: f_refrac=0       ! geometry: 0=Braggref, 1=Laueref
                                              !        2=Braggtrans, 3=Lauetrans
    character(len=sklen)  :: file_refl        ! file created by preprocessor

    ! this stores the data in the preprocessor file
    type(crystalData)     :: crystalData ! type to load xtal info from file_refl

    ! this is defined in module: elasticity, stores the inputs and 
    ! calculated values
    type(crystalElasticity)  :: crystalElasticity ! type to load elastic info

    ! this are inputs specific to bent crystals
    Real(kind=skr)   :: bent_rm=0d0        !meridional bending radius Rm
    Real(kind=skr)   :: bent_rs=0d0        !Sagital bending radius Rs [cm]

    ! these are calculated variables (from crystal_set). Bragg angles
    Real(kind=skr)   :: GRAZE=0d0   !Bragg angle
    Real(kind=skr)   :: THETA_B=0d0 !corrected incidence angle, 
    Real(kind=skr)   :: THETA_B_H=0d0   !corrected reflected angle
    Real(kind=skr)   :: THETA_B_SYM=0d0 !corrected Bragg angle for sym crystal
    Real(kind=skr)   :: THETA_B_COR=0d0 !corrected Bragg angle

    Real(kind=skr)   :: L_EXT_S=0d0     ! Pendellosung length s-pol
    Real(kind=skr)   :: L_EXT_P=0d0     ! idem p-pol
    Real(kind=skr)   :: ASS_FAC=0d0     ! asymmetry factor
    Real(kind=skr)   :: SSR=0d0         ! half-Darwin width s-pol
    Real(kind=skr)   :: SPR=0d0         ! idem p-plo
    Real(kind=skr)   :: Q_MOS=0d0       ! mosaics Q value
    ! bent crystal calculations (e.g., bent_rs will contain the anticlastic)
    Real(kind=skr)   :: bent_m1=0d0        !meridional bending moment
    Real(kind=skr)   :: bent_m2=0d0        !sagittal bending moment
    Real(kind=skr)   :: bent_rm_calc=0d0   !meridional bending radius Rm
    Real(kind=skr)   :: bent_rs_calc=0d0   !Sagital bending radius Rs [cm]
    Real(kind=skr)   :: ml_c=0d0           ! c(sigma) parameter of the ML 
    !Real(kind=skr)   :: pp_beta=0d0        ! beta parameters iof the PP
    !Real(kind=skr)   :: pp_betac=0d0       ! beta-critical parameters iof the PP
    Real(kind=skr)   :: pp_G=0d0           ! G parameter of the PP

    ! these are the 2pi/lambda and photon energy at which the crystal has 
    ! been set (by crystal_set)
    Real(kind=skr)   :: Q_PHOT_SET=0d0 ,PHOT_SET=0d0 

    !vectors
    real(kind=skr),dimension(3)   :: VNOR ! normal to the surface (pointing out)
    real(kind=skr),dimension(3)   :: VPAR ! on the crystal, in diffraction plane
    real(kind=skr),dimension(3)   :: AXIS ! on the crystal, perp to diff plane
    real(kind=skr),dimension(3)   :: BH   ! crystal reciprocal lattice vector
    real(kind=skr),dimension(3)   :: BH0  ! normalized crystal reciprocal lattice vector

    real(kind=skr),dimension(3)   :: VIN_BRAGG ! Bragg direction 
                                               ! (corrected for refraction)
    real(kind=skr),dimension(3)   :: VIN_BRAGG_UNCORR ! Bragg direction 
                                               ! (not corrected for refraction)
    real(kind=skr),dimension(3)   :: VIN_BRAGG_ENERGY
    ! and their respective output vectors
    real(kind=skr),dimension(3)   :: VOUT_BRAGG
    real(kind=skr),dimension(3)   :: VOUT_BRAGG_UNCORR
    real(kind=skr),dimension(3)   :: VOUT_BRAGG_ENERGY
end type oeSetup


!---- Everything is private unless explicitly made public ----!
private 

!---- List of public routines ----!
public :: crystal, scat_direction
public :: crystal_loadCrystalData, plotgle, plotidl, crystal_set_vectors

!---- List of private routines ----!
private :: crystal_printCrystalData
private :: CRYSTAL_FH, CRYSTAL_INFO, CRYSTAL_SET
private :: CRYSTAL_PERFECT, CRYSTAL_MOSAIC,CRYSTAL_BENT_PP,CRYSTAL_BENT_ML
private :: CRYSTAL_BENT_MOMENTS,CRYSTAL_BENT_PP_CALC_BETA,CRYSTAL_BENT_ML_CALC_C

Contains
!
!---- Routines ----!
!
!C
!C       subroutine scat_direction
!C       calculates the output scattering direction applying the
!C       equations k_out_par = k_in_par + BH_par and |k_in| = |k_out|
!C
!C       kk% = ( 2 pi / lambda )  * vv%
!C       vv% direction versor
!C       q_phot = 2 pi / lambda [cm**-1]
!C       bh = ( 2 pi / d_spacing )  normal_to_bragg_planes
!C       vnor = surface_normal
!C       
subroutine scat_direction (vvin1,bh,vnor1,q_phot,vvout)

implicit        none

real(kind=skr),dimension(3),intent(in)   :: vvin1,bh,vnor1
real(kind=skr),intent(in)                :: q_phot
real(kind=skr),dimension(3),intent(out)  :: vvout

real(kind=skr)   :: mod2_kkout_par,mod_kkout_perp,tmp1,tmp2
real(kind=skr),dimension(3)   :: vnor,vvin,kkin,kkin_par,kkin_perp
real(kind=skr),dimension(3)   :: kkout,kkout_par,kkout_perp
real(kind=skr),dimension(3)   :: bh_par,bh_perp
real(kind=skr),dimension(3)   :: vvout1,vvout2,vvoutEwald

vnor = vnor1
vvin = vvin1
Call norm (vnor,vnor)

Call norm (vvin,vvin)
Call scalar (vvin,q_phot,kkin)
Call proj (kkin,vnor,kkin_perp)
Call vector (kkin_perp,kkin,kkin_par)
Call proj (bh,vnor,bh_perp)
Call vector (bh_perp,bh,bh_par)

Call vsum (kkin_par,bh_par,kkout_par)
Call dot(kkout_par,kkout_par,mod2_kkout_par )
mod_kkout_perp = dsqrt(q_phot**2 - mod2_kkout_par)
Call scalar (vnor,mod_kkout_perp,kkout_perp)
! plus sign in  mod_kkout_perp
Call vsum (kkout_par,kkout_perp,kkout)
Call norm (kkout,vvout1)
! minus sign in  mod_kkout_perp
Call vsum (kkout_par,-kkout_perp,kkout)
Call norm (kkout,vvout2)
! select vvout1 or vvout2 
vvoutEwald = bh+vvin*q_phot
call norm(vvoutEwald,vvoutEwald)

call dot(vvout1,vvoutEwald,tmp1)
call dot(vvout2,vvoutEwald,tmp2)
tmp1 = 1.0 - abs(tmp1)
tmp2 = 1.0 - abs(tmp2)
vvout = vvout1
if (tmp2.lt.tmp1) vvout = vvout2

return
end subroutine scat_direction


!C+++
!C       SUBROUTINE    CRYSTAL_FH
!C
!C       PURPOSE         Computes the structute factor of a crystal
!C                       from data stored in crystalData type
!C
!C       ALGORITHM        
!C
!C       MODIFIED        Created by M. Sanchez del Rio  (Feb. 1996)
!C
!C                       Modified by MSR 96/12/19 to include the number
!C                       of coefficients in the f0 parametrization. This
!C                       number is 9 for the CromerMann data and 11 for
!C                       the Waasmaier Kirfel data. Note that the crystal
!C                       file has changed!
!C
!C---
!C 
!C INPUT PARAMETERS:
!C   fUnit: fileUnit
!C     fUnit = [-Inf,3) Not used
!C     fUnit >= 3       Writes info in file unit=fUnit
!C   oe1: oeSetup type with crystal data
!C   PHOT [real]: Photon energy in eV 
!C   THETA [real]: Scattering grazing angle in rads 
!C
!C OUTPUT PARAMETERS:
!C   FH           |
!C   FH_BAR       |
!C   F_0          |    Complex with the returned structure factors 
!C   PSI_H        |    F, Psi and the refraction index.
!C   PSI_HBAR     |
!C   PSI_0        |
!C   REFRAC       |
!C
!C
SUBROUTINE CRYSTAL_FH(fUnit, oe1, PHOT, THETA, &                      ! inputs
                    FH,FH_BAR,F_0,PSI_H,PSI_HBAR,PSI_0,REFRAC,STRUCT) ! outputs 

implicit none

integer(kind=ski),intent(in)::  fUnit
type(oeSetup),intent(in)    ::  oe1
real(kind=skr),intent(in)   ::  PHOT,THETA 
complex(kind=skx),intent(out) :: FH,FH_BAR,F_0
complex(kind=skx),intent(out) :: psi_h,psi_hbar,psi_0,REFRAC,STRUCT

integer(kind=ski):: if0center
integer(kind=ski):: i,j,nener,ierr

real(kind=skr),dimension(NMAXATOMS) :: F0,F1,F2,F2compton

real(kind=skr)   :: r_lam0,ratio,version
real(kind=skr)   :: delta_ref,absorp,absorpNoCompton
complex(kind=skx):: psi_conj,F_0compton
complex(kind=skx):: ci,F(NMAXATOMS)

integer(kind=ski):: i_debug=0

CI = (0.0D0,1.0D0)

if (i_debug.GT.0) then
  write(i_debug,*) '<><> CRYSTAL_FH: fUnit=',fUnit
  write(i_debug,*) '<><> CRYSTAL_FH: oe1%crystalData%d_spacing=',oe1%crystalData%d_spacing
  write(i_debug,*) '<><> CRYSTAL_FH: oe1%crystalData%rn=',oe1%crystalData%rn
endif

!C
!C Computes structure factor at given wavelength and angle.
!C
if (i_debug.EQ.1) then
  write (i_debug,*) '<><>CRYSTAL_FH: working at energy: ',phot
  write (i_debug,*) '<><>CRYSTAL_FH: working at angle [rads] : ',theta
  write (i_debug,*) '<><>CRYSTAL_FH: working at angle [deg] : ',theta*todeg
endif

IF (PHOT.LT.oe1%crystalData%ENERGY(1).OR.PHOT.GT.oe1%crystalData%ENERGY(oe1%crystalData%npoint)) THEN
 write (*,*) 'CRYSTAL_FH: Error: Incoming photon energy is out of range: ',PHOT
 RETURN
END IF

!C
!C Build the fo scattering form factor from its coefficients
!C
ratio = sin(theta)/(toangs/phot)
if (ratio.GT.2) then 
 write (*,*) 'CRYSTAL_FH: ratio sin(theta)/lambda > 2'
 write (*,*) '    ratio : ',ratio
 write (*,*) '    Paramatrizatiog for Fo may fail.'
end if
if (i_debug.gt.0) write (i_debug,*) 'CRYSTAL_FH: ratio is : ',ratio
if0center = (1+oe1%crystalData%nf0coeff)/2

do j=1,oe1%crystalData%nbatom
  f0 (j) = oe1%crystalData%f0coeff(j,if0center)
  do i=1,if0center-1
    f0(j) = f0(j) + oe1%crystalData%f0coeff(j,i) * &
            dexp(-1.0d0*oe1%crystalData%f0coeff(j,i+if0center)*ratio**2)
    if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_FH:  f0 i,j = ',f0(j),i,j
  end do
  if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_FH: F0(j) = ',F0(j)
end do
!C
!C Interpolate for the atomic scattering factor.
!C
DO 299 I = 1, oe1%crystalData%npoint
299         IF (oe1%crystalData%ENERGY(I).GT.PHOT) GO TO 101
!C
I = oe1%crystalData%npoint
101       NENER = I - 1
do j=1,oe1%crystalData%nbatom
 F1(j)= oe1%crystalData%FP(j,NENER) + (oe1%crystalData%FP(j,NENER+1) - & 
      oe1%crystalData%FP(j,NENER)) * (PHOT - oe1%crystalData%ENERGY(NENER)) / & 
      (oe1%crystalData%ENERGY(NENER+1) - oe1%crystalData%ENERGY(NENER))
 F2(j)= oe1%crystalData%FPP(j,NENER) + (oe1%crystalData%FPP(j,NENER+1) - & 
      oe1%crystalData%FPP(j,NENER)) * (PHOT - oe1%crystalData%ENERGY(NENER)) / &
      (oe1%crystalData%ENERGY(NENER+1) - oe1%crystalData%ENERGY(NENER))
 F2compton(j)= &
      oe1%crystalData%FPP(j,NENER)*oe1%crystalData%Fcompton(j,NENER) +&
     (oe1%crystalData%FPP(j,NENER+1)*oe1%crystalData%Fcompton(j,NENER+1) - & 
      oe1%crystalData%FPP(j,NENER)*oe1%crystalData%Fcompton(j,NENER)) * &
      (PHOT - oe1%crystalData%ENERGY(NENER)) / &
      (oe1%crystalData%ENERGY(NENER+1) - oe1%crystalData%ENERGY(NENER))
 if (i_debug.gt.0) then 
   write (i_debug,*) '<><>CRYSTAL_FH: F1 = ',F1(j)
   write (i_debug,*) '<><>CRYSTAL_FH: F2 = ',F2(j)
   write (i_debug,*) '<><>CRYSTAL_FH: Fcompton(j,NENER) = ', &
               oe1%crystalData%Fcompton(j,NENER)
 endif
end do
R_LAM0 = TOCM/PHOT
if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_FH:  R_LAM0  =', R_LAM0

do j=1,oe1%crystalData%nbatom
 F(j)= F0(j) + F1(j) + CI*F2(j)
 if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_FH: F = ',F(j)
end do

!C       
!C FH and FH_BAR are the structure factors for (h,k,l) and (-h,-k,-l).
!C
F_0        = (0.0D0, 0.0D0)
F_0compton = (0.0D0, 0.0D0)
FH         = (0.0D0, 0.0D0)
FH_BAR     = (0.0D0, 0.0D0)
do i=1,oe1%crystalData%nbatom
  FH = FH + (oe1%crystalData%G(i) * F(i) * oe1%crystalData%TEMPER(i) * &
       oe1%crystalData%FRACT(i)) 
  !C From Brennan's private communication:
  !C Changing the temperature doesn't change F0, but it does change
  !C Chi0, due to changing the size of the unit cell, which changes
  !C gamma.
  !C MSR 96/02/14
  
  !
  ! todo: check whether the compton correction affects also f1. 
  !       Now, it only affects the imaginary part f2 but not for the structure
  !       factor, only for the refraction index (and attenuation coeff). 
  F_0 = F_0 + oe1%crystalData%G_0(i) * (oe1%crystalData%atnum(i) + F1(i) + &
        CI * F2(i)) * oe1%crystalData%FRACT(i) 
  F_0compton = F_0compton +  &
               oe1%crystalData%G_0(i) * (oe1%crystalData%atnum(i) + F1(i) + &
               CI * F2compton(i)) * oe1%crystalData%FRACT(i) 
  FH_BAR = FH_BAR + (oe1%crystalData%G_BAR(i) * F(i) * &
           oe1%crystalData%TEMPER(i) * oe1%crystalData%FRACT(i))
end do

STRUCT = sqrt(FH * FH_BAR) 

if (i_debug.gt.0) then 
 write (i_debug,*) '<><>CRYSTAL_FH: FH = ',FH
 write (i_debug,*) '<><>CRYSTAL_FH: FH_BAR = ',FH_BAR
 write (i_debug,*) '<><>CRYSTAL_FH: f_0 = ',f_0
 write (i_debug,*) '<><>CRYSTAL_FH: STRUCT = ',STRUCT
 write (i_debug,*) ' '
 write (i_debug,*) '<><>CRYSTAL_FH: oe1%crystalData%RN = ',oe1%crystalData%RN
 write (i_debug,*) '<><>CRYSTAL_FH: r_lam0 = ',r_lam0
endif

!C
!C   PSI_CONJ = F*( note: PSI_HBAR is PSI at -H position and is
!C   proportional to fh_bar but PSI_CONJ is complex conjugate os PSI_H) 
!C
! srio@esrf.eu 20140116  added sign (see Zachariasen pag 114, eq. 3.95)
psi_h    = -oe1%crystalData%rn*r_lam0**2/pi*fh
psi_hbar = -oe1%crystalData%rn*r_lam0**2/pi*fh_bar
psi_0    = -oe1%crystalData%rn*r_lam0**2/pi*f_0
psi_conj = -oe1%crystalData%rn*r_lam0**2/pi*conjg(fh)

if (i_debug.gt.0) then
  write (i_debug,*) '<><>CRYSTAL_FH: PSI_H = ',PSI_H
  write (i_debug,*) '<><>CRYSTAL_FH: PSI_HBAR = ',PSI_HBAR
  write (i_debug,*) '<><>CRYSTAL_FH: PSI_0 = ',PSI_0
endif

!C
!C computes refractive index.
!C ([3.171] of Zachariasen's book)
!C

!
! values without compton correction
!
REFRAC    = (1.0D0,0.0D0) - R_LAM0**2*oe1%crystalData%RN*F_0/TWOPI
DELTA_REF = 1.0D0 - DREAL(REFRAC)
ABSORP    = 2.0d0 * TWOPI *(-DIMAG(REFRAC)) / R_LAM0
if (i_debug.gt.0) then
  write (i_debug,*) '<><>CRYSTAL_FH: REFRAC NO_COMPTON= ', REFRAC
  write (i_debug,*) '<><>CRYSTAL_FH: DELTA_REF NO_COMPTON= ', DELTA_REF
  write (i_debug,*) '<><>CRYSTAL_FH: ABSORP NO_COMPTON= ', ABSORP
endif

!
! values WITH compton correction
!
absorpNoCompton = absorp
REFRAC = (1.0D0,0.0D0) - R_LAM0**2*oe1%crystalData%RN*F_0compton/TWOPI
DELTA_REF  = 1.0D0 - DREAL(REFRAC)
ABSORP        = 2.0d0 * TWOPI *(-DIMAG(REFRAC)) / R_LAM0
if (i_debug.EQ.1) then
  write (i_debug,*) '<><>CRYSTAL_FH: REFRAC = ', REFRAC
  write (i_debug,*) '<><>CRYSTAL_FH: DELTA_REF = ', DELTA_REF
  write (i_debug,*) '<><>CRYSTAL_FH: ABSORP = ', ABSORP
endif

!C
!C if fUnit > 3, write text to unit fUnit
!C
if (fUnit.gt.3) then 
write(fUnit,*) ' '
write(fUnit,*) '******************************************************'
write(fUnit,*) '           **  DIFF_PAT  v1.8  (24 Mar 2014) **       '
write(fUnit,*) 'Crystal data from file '//trim(oe1%crystalData%filename)
WRITE(fUnit,*) '       at energy    = ',PHOT,' eV'
WRITE(fUnit,*) '                    = ',R_LAM0*1E8,' Angstroms'
WRITE(fUnit,*) '       and at angle = ',THETA*TODEG,' degrees'
WRITE(fUnit,*) '                    = ',THETA,' rads'
write(fUnit,*) '******************************************************'
WRITE(fUnit,*) ' '
do j=1,oe1%crystalData%nbatom
  WRITE(fUnit,*) 'For atom ',j,':'
  WRITE(fUnit,*) '       fo + fp+ i fpp = '
  WRITE(fUnit,*) '        ',f0(j),'+',f1(j),'+ i',f2(j),'= '
  WRITE(fUnit,*) '        ',f0(j)+f1(j)+ci*f2(j)
  WRITE(fUnit,*) '       Z = ',oe1%crystalData%atnum(j)
end do
WRITE(fUnit,*) 'Structure factor F(0,0,0) = ',F_0
WRITE(fUnit,*) 'Structure factor FH = ',FH
WRITE(fUnit,*) 'Structure factor FH_BAR = ',FH_BAR
WRITE(fUnit,*) 'Structure factor F(h,k,l) = sqrt(FH*FH_BAR) ',STRUCT
WRITE(fUnit,*) 'Psi_0  = ',psi_0
WRITE(fUnit,*) 'Psi_H  = ',psi_h
WRITE(fUnit,*) 'Psi_HBar  = ',psi_hbar
WRITE(fUnit,*) '|Psi_H/Psi_HBar|  = ',cdabs(psi_h/psi_hbar)
WRITE(fUnit,*) 'sqrt(Psi_H*Psi_HBar)  = ',sqrt(psi_h*psi_hbar)
WRITE(fUnit,*) '|Psi_H|  = ',cdabs(psi_h)
WRITE(fUnit,*) 'Refraction index = 1 - delta - i*beta :'
WRITE(fUnit,*) '           delta = ',DELTA_REF
WRITE(fUnit,*) '            beta = ',-1.0D0*DIMAG(REFRAC)
WRITE(fUnit,*) 'Absorption coeff = ',ABSORP,' cm^-1'
WRITE(fUnit,*) 'Absorption correction (mu_total/mu_photoelec) = ',ABSORP/absorpNoCompton
WRITE(fUnit,*) ' '
WRITE(fUnit,*) 'e^2/(mc^2)/V = ',oe1%crystalData%rn,' cm^-2'
WRITE(fUnit,*) 'd-spacing = ',oe1%crystalData%d_spacing*1.0E8,' Angstroms'
WRITE(fUnit,*) 'SIN(theta)/Lambda = ',Ratio
WRITE(fUnit,*) ' '
endif
END SUBROUTINE CRYSTAL_FH



!C+++
!C       SUBROUTINE      crystal_loadCrystalData
!C
!C       PURPOSE         Stores the crystallographic data from externam file 
!C                       created with a preprocessor (xop/xcrystal or bragg)
!C
!C       ALGORITHM        
!C
!C       MODIFIED        Created by M. Sanchez del Rio 
!C
!C                       Modified by MSR 96/12/19 to include the number
!C                       of coefficients in the f0 parametrization. This
!C                       number is 9 for the CromerMann data and 11 for
!C                       the Waasmaier Kirfel data. Note that the crystal
!C                       file has changed!
!C                       Modified by MSR 2011/08/16 to port to g95 and 
!C                       integrate in Shadow3
!C
!C---
!C
!C INPUT AND OUTPUT PARAMETERS:
!C   FILEIN : input, the name of the file with data
!C   xtal:    output, the type containing the data loaded from file. 
!C
!C
SUBROUTINE crystal_loadCrystalData (FILEIN,xtal) 

implicit none

character(len=sklen),intent(in) :: FILEIN
type(crystalData),intent(out)   :: xtal


integer(kind=ski)   :: if0center
integer(kind=ski)   :: i,j,nener,ierr
integer(kind=ski)   :: i_debug=0
character(len=sklen):: text


if (i_debug.gt.0) write(i_debug,*) '<><> crystal_loadCrystalData: file is: '//trim(FILEIN)
  OPEN (25,FILE=FILEIN,STATUS='OLD', FORM='FORMATTED',ERR=77)
  read(25,'(A)',err=79)  text
  read(25,*,err=79)  xtal%version, xtal%itype
  if (xtal%itype.NE.1) then 
    write (*,*) 'CRYSTAL_LOADCRYSTALDATA: Error: Data file type not yet implemented.'
   stop
  end if
  read(25,'(A)',err=79)  text
  read(25,*,err=79)  xtal%RN,xtal%d_spacing  !,TEMPER
  read(25,'(A)',err=79)  text
  read(25,*,err=79)  xtal%nbatom
  if (xtal%nbatom.GT.NMAXATOMS) then
    write (*,*) 'CRYSTAL_LOADCRYSTALDATA: Error: Maximum number of atoms allowad: ',NMAXATOMS
  end if
  
  read(25,'(A)',err=79) text
  read(25,*,err=79)  (xtal%ATNUM(I),i=1,xtal%nbatom)
  
  read(25,'(A)',err=79) text
  read(25,*,err=79)  (xtal%FRACT(I),i=1,xtal%nbatom)
  
  read(25,'(A)',err=79) text
  read(25,*,err=79)  (xtal%TEMPER(I),i=1,xtal%nbatom)
  
  read(25,'(A)',err=79) text
  read(25,*,err=79) (xtal%G_0(i),i=1,xtal%nbatom)
  
  read(25,'(A)',err=79)  text
  do i=1,xtal%nbatom
    read(25,*,err=79)  xtal%G(I)
    read(25,*,err=79)  xtal%G_BAR(I)
  end do
  
  read(25,'(A)',err=79) text
  do i=1,xtal%nbatom
    read(25,*,err=79)  xtal%nf0coeff,(xtal%f0coeff(i,j),j=1,xtal%nf0coeff)
  end do
  read(25,'(A)',err=79)  text
  read(25,*,err=79)  xtal%NPOINT
  read(25,'(A)',err=79)  text
  do I = 1, xtal%NPOINT
    read(25,*,err=79)  xtal%energy(i)
    if ( (xtal%version-2.39).GT.0) then
        ! new version (since 2.4) implements compton correction factor
        ! mu = mu_photoelectric*Fcompton
        ! where Fcompton=mu_total/mu_photoelectric
        ! and mu_total is usually mu_photoelectric+mu_compton
        ! (elastic may be also added)
        do j=1,xtal%nbatom
           read (25,*,err=79) xtal%FP(j,i),xtal%FPP(j,i),xtal%Fcompton(j,i)
        end do
    else
        do j=1,xtal%nbatom
           read (25,*,err=79) xtal%FP(j,i),xtal%FPP(j,i)
           xtal%Fcompton(j,i) = 1.0d0
        end do
    endif
  END DO 
  xtal%filename = FILEIN

CLOSE (25)
RETURN

77 CONTINUE
79 CONTINUE
   print *,'crystal_loadCrystalData: Error reading file: '//trim(FILEIN)
   stop
END SUBROUTINE crystal_loadCrystalData


!C+++
!C       SUBROUTINE    crystal_printCrystalData
!C
!C       PURPOSE         Print the stored crystallographic data from
!C                       preprocessor (xop/xcrystal or bragg)
!C
!C       ALGORITHM        
!C
!C       MODIFIED        Created by M. Sanchez del Rio 
!C
!C
!C---
SUBROUTINE crystal_printCrystalData (xtal) 

implicit none

type(crystalData),intent(in) :: xtal
integer(kind=ski)            :: i,j

  print *,' '
  print *,'============== in crystal_printCrystalData =================='
  print *,'NMAXATOMS=',NMAXATOMS
  print *,'NMAXENER=',NMAXENER
  print *,'RN=',xtal%rn
  print *,'d_spacing=',xtal%d_spacing
  print *,'ATNUM=',xtal%atnum(1:xtal%NBATOM)
  print *,'TEMPER=',xtal%temper(1:xtal%NBATOM)
  print *,'NBATOM,=',xtal%nbatom
  print *,'G=',xtal%g(1:xtal%NBATOM)
  print *,'G_BAR=',xtal%g_bar(1:xtal%NBATOM)
  print *,'G_0=',xtal%g_0(1:xtal%NBATOM)
  print *,'NF0COEFF=',xtal%nf0coeff
  print *,'F0COEFF=',xtal%f0coeff(1:xtal%NBATOM,1:11)
  print *,'NPOINT=',xtal%npoint
  print *,'FRACT=',xtal%fract(1:xtal%NBATOM)
  print *,'filename=',trim(xtal%filename)
  print *,"energy  atom_index  FP  FPP  "
  print *,'============================================================='
  print *,' '

  DO I = 1, xtal%NPOINT
    do j=1,xtal%nbatom
       print *,xtal%energy(i),j, xtal%FP(j,i),xtal%FPP(j,i)
    end do
  END DO 
  RETURN

END SUBROUTINE crystal_printCrystalData

!C+++
!C       SUBROUTINE    CRYSTAL
!C
!C       PURPOSE       Computes the reflectivity of a symmetric Bragg crystal 
!C                     according to the dynamic theory of x-ray diffraction.
!C
!C       ALGORITHM     Reference B.E.Warren, X-Ray Diffraction, 
!C                     Addison-Wesley 1969.  See also M.J.Bedzyk, G.Materlik 
!C                     and M.V.Kovalchuk, Phy. Rev. B30, 2453(1984).
!C                     For mosaic crystal reflectivity see Zachariasen, 
!C                     Theory of x-ray diffraction in crystals, Dover (1966)
!C                     formula 4.24, and Bacon and Lowde, Acta Crystall.1
!C                     pag 303 (1948) formula 17. 
!C
!C       MODIFIED      July 1989, M. Sanchez del Rio for asymmetry part,
!C                     July 1990, mosaic part.
!C                     August 1992, laue part.
!C                     Dec 2012, cleaned
!C
!C---
!C
!C kwhat: flag:
!C       <0 read crystal file (call crystal_loadCrystalData) (note that 
!C                             thio option changes oe1)
!C       0 set mode (call crystal_set)
!C       1 info mode (call crystal_info)
!C       2 reflectivity mode (call crystal_{perfect,mosaic})

SUBROUTINE CRYSTAL (KWHAT, &               ! input
           oe1, &                          ! input (also output if kwhat<1)
           Q_PHOT, &                       ! inputs (kwhat>=0)
           VIN, VOUT, BH, SURFNOR, &       ! inputs (kwhat=2)
           R_S, R_P,PHASE_S, PHASE_P, L_EXT_S, L_EXT_P) !, & !  output(kwhat=2)
           !THETA_B, SSR, SPR)                  !  output(kwhat=2)

implicit none
                               
integer(kind=ski),intent(in) ::  KWHAT
type(oeSetup),intent(inout)  ::  oe1
real(kind=skr),intent(in)    ::  Q_PHOT                              !Arguments
real(kind=skr),dimension(3),intent(in) ::  VIN, VOUT, BH, SURFNOR
real(kind=skr),intent(out)   ::  R_S, R_P,PHASE_S, PHASE_P
real(kind=skr),intent(out)   ::  L_EXT_S, L_EXT_P

real(kind=skr)   ::  ass_fac
real(kind=skr)   ::  phot,r_lam0
real(kind=skr)   ::  graze
real(kind=skr)   ::  absorp    
real(kind=skr)   ::  q_mos
real(kind=skr)   :: theta_b_h,theta_b_cor,theta_b_sym
complex(kind=skx):: F_0,FH,FH_BAR,REFRAC,STRUCT
complex(kind=skx):: psi_h,psi_hbar,psi_0

integer(kind=ski)::  i_debug=0


! first the usual mode (kwhat=2) to be fast

if (kwhat.eq.2) then 
    !C***************************************************
    !C If flag is = 2, gets reflectivity
    !C***************************************************
    
    PHOT = Q_PHOT/TWOPI*TOCM
    R_LAM0 = TWOPI/Q_PHOT
    !THETA_B = oe1%THETA_B

    ! these are needed in output 
    L_EXT_S = oe1%L_EXT_S
    L_EXT_P = oe1%L_EXT_P

    GRAZE = ASIN(R_LAM0/oe1%crystalData%d_spacing/2.0D0)
    
    if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL: working at energy: ',phot
    !!C
    !!C call crystal_fh to get the structure factor
    !!C
    if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL: calling crystal_fh: ',kwhat
    Call crystal_fh (KWHAT, oe1, PHOT, GRAZE , &
           FH,FH_BAR,F_0,PSI_H,PSI_HBAR,PSI_0,REFRAC,STRUCT)
    
    if (i_debug.gt.0) then
       write(i_debug,*) '<><>CRYSTAL: back crystal_fh, oe1%crystalData%rn: ',oe1%crystalData%rn
       write(i_debug,*) 'CRYSTAL: THETA_B=',oe1%THETA_B
       write(i_debug,*) 'CRYSTAL: THETA_B[deg]=',oe1%THETA_B*TODEG
       write(i_debug,*) 'CRYSTAL: L_EXT_S=',oe1%L_EXT_S
    end if
    
    !C
    !C Main calculation
    !C
      IF (oe1%F_MOSAIC.EQ.1) THEN   ! mosaic crystal
       Q_MOS = pi**2*abs(psi_h*psi_hbar)/r_lam0/sin(2*graze)
       !note that crystal)mosaic overwrite the extinction with 
       !secondary extinction
       CALL CRYSTAL_MOSAIC   ( Q_PHOT, VIN, BH, SURFNOR, oe1, &
          Q_MOS,REFRAC, R_S, R_P,PHASE_S, PHASE_P, L_EXT_S, L_EXT_P)
      ELSEif(oe1%F_MOSAIC.EQ.0)then ! perfect crystal
       call CRYSTAL_PERFECT         (Q_PHOT, VIN, VOUT, BH, SURFNOR,oe1, &
          PSI_0,PSI_H,PSI_HBAR, R_S, R_P,PHASE_S, PHASE_P)
      ELSEif(oe1%F_MOSAIC.EQ.2)then ! bent, multilamellar
       call CRYSTAL_BENT_ML (Q_PHOT, VIN, VOUT, BH, SURFNOR, oe1, &
          PSI_0,PSI_H,PSI_HBAR,REFRAC, STRUCT, FH,FH_BAR,R_S, R_P,PHASE_S, PHASE_P)
      Elseif(oe1%F_MOSAIC.EQ.3)then ! bent, penning-polder
        call CRYSTAL_BENT_PP (Q_PHOT, VIN, VOUT,BH, SURFNOR, oe1, &
        PSI_0,PSI_H,PSI_HBAR,STRUCT, REFRAC,R_S, R_P,PHASE_S, PHASE_P)
      Elseif(oe1%F_MOSAIC.EQ.4)then
      !call bent tt
      END IF
    
return
end if

! then the unusual mode (kwhat < 2) 

select case(kwhat)
  case(0)  
         CALL CRYSTAL_SET (Q_PHOT,oe1) 
  case(1)  
         CALL CRYSTAL_INFO (oe1) 
  case default   ! load file
         Call crystal_loadCrystalData (oe1%FILE_REFL,oe1%crystalData)
         if (i_debug.gt.0) Call crystal_printCrystalData (oe1%crystalData)
end select

return
 
END SUBROUTINE CRYSTAL

!C+++
!C       SUBROUTINE       CRYSTAL_INFO
!C
!C       PURPOSE          Prints (on the screen and on fortran unit 23) some 
!C                        info from the crystal set at a given photon energy
!C
!C       ALGORITHM        
!C
!C       MODIFIED         M. Sanchez del Rio
!C
!C---

SUBROUTINE CRYSTAL_INFO (oe1) 

implicit none
                               
type(oeSetup),intent(inout)    :: oe1

real(kind=skr)   :: ssr
real(kind=skr)   :: spr
real(kind=skr)   :: THETA_B
integer(kind=ski):: i
real(kind=skr)   :: q_phot,phot,r_lam0
real(kind=skr)   :: graze,sin_gra,ass_fac
real(kind=skr)   :: absorp   
real(kind=skr)   :: q_mos
real(kind=skr)   :: tmax_mos
real(kind=skr)   :: texts,textp
complex(kind=skx):: F_0,FH,FH_BAR,REFRAC,STRUCT
complex(kind=skx):: psi_h,psi_hbar,psi_0
real(kind=skr)   :: theta_b_sym,theta_b_h,theta_b_cor
real(kind=skr)   :: sigma_gamma,biga_mosaic,omega,kpolsqrt

integer(kind=ski):: i_debug=0

q_phot = oe1%q_phot_set
phot = oe1%phot_set
!C
!C Computes reflectivities at given wavelength and angle.
!C

!!  PHOT        = Q_PHOT/TWOPI*TOCM
!!  R_LAM0         = TWOPI/Q_PHOT
!!  SIN_GRA = R_LAM0/oe1%crystalData%d_spacing/2.0D0
!!  GRAZE = ASIN(SIN_GRA)
GRAZE = oe1%graze

!todo: be sure that crystal_set has been called before calling crystal_info
!CALL CRYSTAL_SET (Q_PHOT,oe1) !, &
!         !THETA_B,THETA_B_H,THETA_B_SYM, THETA_B_COR, L_EXT_S,L_EXT_P,ASS_FAC,SSR,SPR,Q_MOS)
THETA_B = oe1%THETA_B
THETA_B_H = oe1%THETA_B_H
THETA_B_SYM = oe1%THETA_B_SYM
THETA_B_COR = oe1%THETA_B_COR
!L_EXT_S = oe1%L_EXT_S
!L_EXT_P = oe1%L_EXT_P
ASS_FAC = oe1%ASS_FAC
SSR = oe1%SSR
SPR = oe1%SPR
Q_MOS = oe1%Q_MOS

if (i_debug.gt.0) then
 write (i_debug,*) '<><>CRYSTAL_INFO ssr: ',ssr
 write (i_debug,*) '<><>CRYSTAL_INFO spr: ',spr
 write (i_debug,*) '<><>CRYSTAL_INFO theta_b: ',theta_b
end if


!C***************************************************
!C display results on screen (unit=6) and in unit 27 
!C***************************************************
!C

!C
!C Calculate mosaic parameters
!C
if (oe1%f_mosaic.EQ.1) then
 if (oe1%f_refrac.eq.0) then
   biga_mosaic = oe1%thickness*absorp/sin(torad*graze)    !bragg
 else if (oe1%f_refrac.eq.1) then
   biga_mosaic = oe1%thickness*absorp/cos(torad*graze)    !laue(alpha=90)
 else 
   write (*,*) 'CRYSTAL_INFO Error: Option not implemented.'
   stop
 endif
 omega   = (1/sqrt(twopi))*(1/oe1%spread_mos)
!C
!C Q_mos is the mosaic q variable (Zachariasen's [4.33b])
!C Sigma_Gamma refers to Zachariasen's formula 4.33a
!C
!q_mos = pi**2*abs(psi_h*psi_hbar)/r_lam0/sin(2*graze)
  sigma_gamma =  omega*q_mos
  kpolsqrt = (cos(2*graze))**2

 tmax_mos = cos(graze)*dexp(1+2*sigma_gamma/absorp)/2/sigma_gamma/absorp/absorp
!C         texts= r_lam0**2/pi/sin(2*graze)/ssr/2/d_spacing
!C         textp= r_lam0**2/pi/sin(2*graze)/spr/2/d_spacing
 texts= oe1%L_EXT_S*sin(GRAZE)
 textp= oe1%L_EXT_P*sin(GRAZE)
end if

DO i=6,23,17       ! write output in units 6(screen) and 23
  WRITE(i,*) ' '
  call crystal_fh (i, oe1, &
       PHOT, GRAZE,FH,FH_BAR,F_0,PSI_H,PSI_HBAR,PSI_0,REFRAC,STRUCT)

  WRITE(i,*) 'Photon energy = ',Q_PHOT/TWOPI*TOCM ,' eV'
  WRITE(i,*) '       wavelength = ',TWOPI/Q_PHOT*1D8,' Angstroms'
  WRITE(i,*) 'Theta (grazing uncorrected Bragg angle)   = ',GRAZE,' rads'
  WRITE(i,*) '                                          = ',GRAZE*TODEG,' degrees'
  WRITE(i,*) 'Symmetric Bragg angle (corrected) = ', THETA_B_SYM,' rad'
  WRITE(i,*) '                                  = ', THETA_B_SYM*TODEG,' degrees'
  WRITE(i,*) 'Angle (H,kin) = Bragg angle (corrected) including asymmetry = ', &
              THETA_B_COR, 'rad'
  WRITE(i,*) '                                                            = ', &
              THETA_B_COR*TODEG,' degrees'
  WRITE(i,*) 'Asymmetry angle alpha= ',oe1%a_bragg*todeg,' degrees'
  WRITE(i,*) 'Asymmetry factor b=  ',ass_fac
  WRITE(i,*) ' '
!C  Following line modified by T. Jach, 10.10.2001
  WRITE(i,*) 'Extinction lengths and depths given here are for amplitude (for intensities, '
  WRITE(i,*) '   consider half value).'
  WRITE(i,*) 'DEPTH values are measured perpendicular to the crystal surface and LENGTH'
  WRITE(i,*)  '   values are along the incident beam path.'
  WRITE(i,*) 'Extinction length (amplitude) sigma= ', oe1%L_EXT_S*1D4/abs(oe1%vin_bragg(3)),' microns'
!C  Following line modified by T. Jach, 10.10.2001
  WRITE(i,*) 'Extinction length (amplitude) pi= ', oe1%L_EXT_P*1D4/abs(oe1%vin_bragg(3)),' microns'
!C  Following line modified by T. Jach, 10.10.2001 
  WRITE(i,*) 'Extinction depth (amplitude) sigma= ',oe1%L_EXT_S*1D4,' microns'
!C  Following line modified by T. Jach, 10.10.2001 
  WRITE(i,*) 'Extinction depth (amplitude) pi= ', oe1%L_EXT_P*1D4,' microns'
!C  Following line modified by  srio 26.10.2001
  WRITE(i,*) 'Pendellosung period (amplitude) sigma = Extinction depth sigma times 2 pi = ', oe1%L_EXT_S*1D4*2D0*PI, ' microns'     
!C  Following line added  srio 26.10.2001
  WRITE(i,*) 'Pendellosung period (amplitude) pi = Extinction depth pi times 2 pi = ', oe1%L_EXT_P*1D4*2.D0*PI, ' microns'     
WRITE(i,*) ' '

 if (oe1%f_bragg_a.eq.0) then

   WRITE(i,*) 'width of diffraction profile  s-pol  =  ', 2.0D0*SSR*1.0D+6,' microradians'
   WRITE(i,*) '                                =  ', 2.0D0*SSR*TODEG*3600D0,' arc sec'
   WRITE(i,*) 'width of diffraction profile  p-pol  =  ', 2.0D0*SPR*1.0D+6,' microradians'
   WRITE(i,*) '                                =  ', 2.0D0*SPR*TODEG*3600D0,' arc sec'
 else
   WRITE(i,*) 'The width of the diffraction profile as a function of incident angle is'
   WRITE(i,*)  'width for s-pol  = ', 2.0D0*SSR*1.0D+6,' microradians'
   WRITE(i,*)  '                     =  ', 2.0D0*SSR*TODEG*3600D0,' arc sec'
   WRITE(i,*)  'width for p-pol  = ', 2.0D0*SPR*1.0D+6,' microradians'
   WRITE(i,*)  '                     =  ', 2.0D0*SPR*TODEG*3600D0,' arc sec'
   !!WRITE(i,*)  'Incident Grazing Angle = ', (oe1%a_bragg+graze)*todeg,' deg'
   WRITE(i,*)  'Incident Corrected Angle = ', (theta_b)*todeg,' deg'
   WRITE(i,*)   ' '
   WRITE(i,*)   'The width of the diffraction profile as a function of reflected angle is'
   WRITE(i,*)  'width for s-pol  = ', 2.0D0*SSR*1.0D+6*(abs(ass_fac)),' microradians'
   WRITE(i,*)  '                     = ', 2.0D0*SSR*(abs(ass_fac))*TODEG*3600D0,' arc sec'
   WRITE(i,*)  'width for p-pol  = ', 2.0D0*SPR*1.0D+6*(abs(ass_fac)),' microradians'
   WRITE(i,*)  '                     = ', 2.0D0*SPR*(abs(ass_fac))*TODEG*3600D0,' arc sec'
   !!WRITE(i,*)  'Reflected Grazing Angle = ', (graze-oe1%a_bragg)*todeg,' deg'
   WRITE(i,*)  'Reflected Corrected Angle = ', (theta_b_h)*todeg,' deg'
 endif
END DO
i = 23 ! only to file
if (oe1%f_mosaic.EQ.1) then
 write (i,*) '  '
 write (i,*) '***********  MOSAIC PARAMETERS  ***************'
 write (i,*) '  '
 write (i,*) 'spread_mos= ',2.35*oe1%spread_mos/TORAD ,' deg fwhm'
 write (i,*) 'true absorp depth = ',sin(graze)/absorp*1d4, ' microns'
 write (i,*) 'true absorp length = ',1d4/absorp,' microns'
 write (i,*) 'peak thickness = ',tmax_mos,' cm'
 write (i,*) ' '
 write (i,*) 'For s-polarization we have: '
 write (i,*) '   Q      = ',q_mos,' cm^-1 '
!C         write (i,*) '   Peak reflectivity = ',r_s
 write (i,*) '   Primary extinction:'
 write (i,*) '      mu =',sin(graze)/texts,'cm^-1'
 write (i,*) '      length =',texts/sin(graze)*1d4,' microns'
 write (i,*) '      depth =',texts*1d4,' microns'
 write (i,*) '   Secondary extinction:'
 write (i,*) '      mu =',sigma_gamma,'cm^-1'
 write (i,*) '      length =',1d4/sigma_gamma,' microns'
 write (i,*) '      depth =',sin(graze)/sigma_gamma*1d4, ' microns'

 write (i,*) ' '
 write (i,*) 'For p-polarization we have: '
 write (i,*) '   Q     = ',q_mos*kpolsqrt,' cm^-1 '
!C         write (i,*) '   Peak reflectivity = ',r_p
 write (i,*) '   Primary extinction:'
 write (i,*) '      mu =',sin(graze)/textp,'cm^-1'
 write (i,*) '      length =',textp/sin(graze)*1d4,' microns'
 write (i,*) '      depth =',textp*1d4,' microns'
 write (i,*) '   Secondary extinction:'
 write (i,*) '      mu =',sigma_gamma*kpolsqrt,'cm^-1'
 write (i,*) '      length =',1d4/sigma_gamma/kpolsqrt, ' microns'
 write (i,*) '      depth =', sin(graze)/sigma_gamma/kpolsqrt*1d4, ' microns'
 write (i,*) '  '
endif

write (i,*) '  '
write (i,*) '  '
write (i,*) '***********  VECTOR AND ANGLES  ***************'
write (i,*) '  '
write (i,*) 'Reciprocal space hkl along BH: ',oe1%BH
write (i,*) 'Reciprocal space hkl along BH (normalized): ',sV(oe1%BH0)
write (i,*) 'Normal to the crystal VNOR:                 ',sV(oe1%VNOR)
write (i,*) '  '
write (i,*) 'Incident directions matching Bragg angle:'
write (i,*) '    VIN_BRAGG_UNCORR (Uncorrected): ',sV(oe1%VIN_BRAGG_UNCORR)
write (i,*) '    VIN_BRAGG          (Corrected): ',sV(oe1%VIN_BRAGG)
write (i,*) '    VIN_BRAGG_ENERGY              : ',sV(oe1%VIN_BRAGG_ENERGY)
write (i,*) 'Reflected directions matching Bragg angle:'
write (i,*) '   VOUT_BRAGG_UNCORR (Uncorrected): ',sV(oe1%VOUT_BRAGG_UNCORR)
write (i,*) '   VOUT_BRAGG          (Corrected): ',sV(oe1%VOUT_BRAGG)
write (i,*) '   VOUT_BRAGG_ENERGY              : ',sV(oe1%VOUT_BRAGG_ENERGY)
write (i,*) 'Gamma0 =  VIN_BRAGG(3) = ',oe1%VIN_BRAGG(3)
write (i,*) 'GammaH = VOUT_BRAGG(3) = ',oe1%VOUT_BRAGG(3)
write (i,*) 'b~Gamma0/GammaH = ',oe1%VIN_BRAGG(3)/oe1%VOUT_BRAGG(3)
write (i,*) '  '
write (i,*) 'angle between VIN_BRAGG_UNCORR and VNOR [deg]:      ', &
            todeg*acos(dot_product(oe1%vin_bragg_uncorr,oe1%vnor))
write (i,*) 'angle between VOUT_BRAGG_UNCORR and VNOR [deg]:     ', &
            todeg*acos(dot_product(oe1%vout_bragg_uncorr,oe1%vnor))

write (i,*) 'angle between VIN_BRAGG and VNOR [deg]:             ', &
            todeg*acos(dot_product(oe1%vin_bragg,oe1%vnor))
write (i,*) 'angle THETA_B: ',todeg*oe1%theta_b
write (i,*) 'angle between VOUT_BRAGG and VNOR [deg]:            ', &
            todeg*acos(dot_product(oe1%vout_bragg,oe1%vnor))
write (i,*) 'angle THETA_B_H [deg] =: ',todeg*oe1%theta_b_h

write (i,*) 'angle between VIN_BRAGG and BH [deg]:               ', &
            todeg*acos(dot_product(oe1%vin_bragg,oe1%bh0))
write (i,*) 'angle THETA_B_COR [deg] =: ',todeg*oe1%theta_b_cor

write (i,*) 'angle between VIN_BRAGG_UNCORR and BH [deg]:        ', &
            todeg*acos(dot_product(oe1%vin_bragg_uncorr,oe1%bh0))
write (i,*) 'angle GRAZE [deg]:    ',todeg*oe1%GRAZE
write (i,*) 'angle alphaX [deg]:   ',todeg*oe1%a_bragg
! alphaX+alpha=180
write (i,*) 'angle alpha [deg]:    ',180d0-todeg*oe1%a_bragg
! alpha=chi+90
write (i,*) 'angle chi [deg]:      ', 90d0-todeg*oe1%a_bragg

if (oe1%f_mosaic.GT.1) then
   write (i,*) '  '
   write (i,*) '  '
   write (i,*) '***********  ELASTICITY PARAMETERS  ***************'
   write (i,*) '  '
   call elasticity_report(oe1%crystalElasticity,i)
   call crystal_bent_moments(oe1,i)
end if

if (oe1%f_mosaic.EQ.2) then
   write (i,*) '  '
   write (i,*) '  '
   write (i,*) '***********  BENT CRYSTAL PARAMETERS (MULTILAMELAR MODEL)  ***************'
   write (i,*) '  '
   call crystal_bent_ml_calc_c(oe1,psi_h,i)
endif

if (oe1%f_mosaic.EQ.3) then
   write (i,*) '  '
   write (i,*) '  '
   write (i,*) '***********  BENT CRYSTAL PARAMETERS (PENNING-POLDER MODEL) ***************'
   write (i,*) '  '
   call crystal_bent_pp_calc_beta(oe1,psi_h,i)
endif

END SUBROUTINE CRYSTAL_INFO

!
! Auxiliary function to print vectors easily
!
FUNCTION sV(vector) !,myformat)
implicit none 
real(kind=skr),dimension(3),intent(in) :: vector
!character(len=sklen),intent(in)        :: myformat
character(len=sklen)                   :: sV
write(sV,111) vector(1),vector(2),vector(3)
return 
111 format('(',f12.8,',',f12.6,',',f12.6,')') 
END FUNCTION sV


!C+++
!C       SUBROUTINE      CRYSTAL_SET
!C
!C       PURPOSE         Calculates useful parameters of the crystal
!C                         (Bragg angle corrected for refraction, widths
!C                          of Darwin profiles, etc.)
!C
!C       ALGORITHM        Reference B.E.Warren, X-Ray Diffraction, 
!C
!C       MODIFIED         M. Sanchez del Rio
!C 
!C                       2012-12-07 changed to fully use Zachariasen equations. 
!C                                  changed "default" laue case to be "onto" 
!C                                  bragg planes
!C                                  outputs are saven in oeSetup structure
!C
!C
!C---

SUBROUTINE CRYSTAL_SET (Q_PHOT, oe1)  

implicit none
                               
real(kind=skr),intent(in)     :: Q_PHOT   !Arguments
type(oeSetup),intent(inout)   :: oe1


real(kind=skr)   :: THETA_B,theta_b_h,theta_b_sym, theta_b_cor
real(kind=skr)   :: L_EXT_S, L_EXT_P
real(kind=skr)   :: ass_fac
real(kind=skr)   :: ssr,spr,q_mos
real(kind=skr)   :: DELTA_REF
real(kind=skr)   :: graze,sin_gra,theta_o,theta_h
real(kind=skr)   :: gammaf, ppol
real(kind=skr)   :: phot,r_lam0
real(kind=skr)   :: absorp,tmp
real(kind=skr)   :: gamma_0,gamma_H
complex(kind=skx):: F_0,FH,FH_BAR,STRUCT,REFRAC
complex(kind=skx):: psi_h,psi_hbar,psi_0

integer(kind=ski):: i,i_debug=0

!C
!C Computes reflectivities at given wavelength and angle.
!C
PHOT    = Q_PHOT/TWOPI*TOCM
R_LAM0  = TWOPI/Q_PHOT
SIN_GRA = R_LAM0/oe1%crystalData%d_spacing/2.0D0
GRAZE   = ASIN(SIN_GRA)

if (i_debug.gt.0) then
  write (i_debug,*) '<><>CRYSTAL_SET energy PHOT=',phot
  write (i_debug,*) '<><>CRYSTAL_SET A_BRAGG [deg] : ',todeg*oe1%a_bragg
  write (i_debug,*) '<><>CRYSTAL_SET GRAZE [deg] : ',todeg*graze
endif

! save uncorrected Bragg angle
oe1%GRAZE       = GRAZE
oe1%Q_PHOT_SET  = Q_PHOT

! now call crystal_set_vectors to define the incident and reflected directions, 
! but only for the uncorrected directions (we do not know yet the corrected 
! angles)

call crystal_set_vectors(oe1)

! calculate gamma (note that our director cosines have oppotite
! sign to those of Zachariasen, as he uses inwards vnor)
gamma_0 = dot_product(oe1%vin_bragg_uncorr,oe1%vnor)
gamma_h = dot_product(oe1%vout_bragg_uncorr,oe1%vnor)
! calculate angles
theta_o = acos(gamma_0)
theta_h = acos(gamma_H)

! ass_fac is b=gamma0/gammaH=vin(3)/vout(3)
!ass_fac = oe1%vin_bragg_uncorr(3)/oe1%vout_bragg_uncorr(3)
ass_fac = gamma_0/gamma_H

!C
!C call crystal_fh to get the structure factor
!C
if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_SET: Calling crystal_fh'
Call crystal_fh (1, oe1, PHOT, GRAZE , &
                 FH,FH_BAR,F_0,PSI_H,PSI_HBAR,PSI_0, &
                 REFRAC,STRUCT)
if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_SET back crystal_fh, oe1%crystalData%rn: ',oe1%crystalData%rn

STRUCT = sqrt(FH * FH_BAR) 
if (i_debug.gt.0) then
  write (i_debug,*) '<><>CRYSTAL_SET FH = ',FH
  write (i_debug,*) '<><>CRYSTAL_SET f_0 = ',f_0
  write (i_debug,*) '<><>CRYSTAL_SET FH_BAR = ', FH_BAR
  write (i_debug,*) '<><>CRYSTAL_SET STRUCT = ', STRUCT
  write (i_debug,*) '<><>CRYSTAL_SET PSI_H = ', PSI_H
  write (i_debug,*)  '<><>CRYSTAL_SET PSI_HBAR = ', PSI_HBAR
  write (i_debug,*) '<><>CRYSTAL_SET PSI_0 = ', PSI_0
endif

DELTA_REF  = 1.0D0 - DREAL(REFRAC)
ABSORP        = 2.0d0 * TWOPI *(-DIMAG(REFRAC)) / R_LAM0

!C
!C Calculates the reflection angles and other useful parameters
!C

! srio@esrf.eu  2012-12-07
! The angles are (see our paper)
! alphaX = alpha as defined in XOP (angle from crystal surface to Bragg plane, >0 if cw)
! alpha = alpha "mathematically correct" (angle from crystal surface (axis y) to Bragg 
! plane in the direction that H points
!
! chi = angle between H and the normal to the surface
! alpha = chi + 90
! alpha + alphaX = 180
! therefore chi = 90-alphaX
!
! angles (in mathematical sense measured from axis onto the crystal surface (i.e., "y"))
! angles are >0 if ccw

If((oe1%f_refrac .eq.1) .or. (oe1%f_refrac .eq.3)) then ! Laue symmetric (b=1)
  THETA_B_SYM=graze
else ! Bragg symmetric (b=-1)
  THETA_B_SYM = graze - Dreal(PSI_0)/sin(2*graze)
end if
if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_SET THETA_B_SYM **new** [deg]: ',THETA_B_SYM*todeg

!
! This is the Bragg corrected angle in the general (asymmetric) case
!
THETA_B_COR = graze + (1d0-ass_fac)/ass_fac*Dreal(PSI_0)/2d0/sin(2d0*graze)
if (i_debug.gt.0) write (i_debug,*) '<><>CRYSTAL_SET THETA_B_COR **new** [deg]: ',THETA_B_COR*todeg

! this is the (unfigned) factor that converts structure factor F into Psi
gammaf = oe1%crystalData%RN*(R_LAM0**2)/PI
! p-polarization factor
PPOL = ABS(COS(2.0D0*GRAZE))

!
! Extinction lengths (along the incident beam)  and 
! depths (along the surface normal)
!

!C L_EXT_S = R_LAM0*cos(theta_b)/gammaf/abs(fh)
!C The formula above is correct for symmetric case (Bragg and Laue). 
!C Bug fixed 14 Set 2001 srio@esrf.fr
!C The extinction length was only for Laue symmetric case.
!C Now, the general formula is implemented. 
!C Thanks to T. Jach who pointed out the problem, and previous
!C remarks by J. Sutter, E. Alp and R. Dejus.
!C
!C Following line modified by T. Jach 10/10/01
!C L_EXT_S = R_LAM0*sin(theta_o)/gammaf/abs(fh)/sqrt(abs(ass_fac))/PI
!C

!TODO: Note that this is the extinction depth in AMPLITUDE. For intensity
!     take 1/2 of it.
!TODO: check cdabs(psi_h) or sqrt(psi_h*psi_hbar)
L_EXT_S = R_LAM0*abs(gamma_0)/cdabs(psi_h)/sqrt(abs(ass_fac))/PI
L_EXT_P = L_EXT_S/PPOL

q_mos = pi**2*abs(psi_h*psi_hbar)/r_lam0/sin(2*graze)

!
! Half Darwin width = SSR 
! in Zachariasen formalism: the half-width angular value corresponding to y=1
!
!ssr = real(sqrt(psi_h*psi_hbar)/sqrt(abs(ass_fac))/sin(2*graze))
ssr = cdabs(psi_h)/sqrt(abs(ass_fac))/sin(2*graze)
spr = ssr*ppol
if (i_debug.gt.0) then
   write (i_debug,*) '<><>CRYSTAL_SET SSR **new** : ',ssr
   write (i_debug,*) '<><>CRYSTAL_SET SPR **new** : ',spr
endif

!
! so far, store results in oeSetup structure
!
oe1%THETA_B_SYM = THETA_B_SYM
oe1%THETA_B_COR = THETA_B_COR
oe1%L_EXT_S     = L_EXT_S
oe1%L_EXT_P     = L_EXT_P
oe1%ASS_FAC     = ASS_FAC
oe1%SSR         = SSR
oe1%SPR         = SPR
oe1%Q_MOS       = Q_MOS
oe1%PHOT_SET    = PHOT

! now calculate again the directions but now also the 
! corrected (for refraction) directions
call crystal_set_vectors(oe1)

theta_b = acos(dot_product(oe1%vin_bragg,oe1%vnor))
theta_b_h = acos(dot_product(oe1%vout_bragg,oe1%vnor))

oe1%THETA_B     = THETA_B
oe1%THETA_B_H   = THETA_B_H


! now calculate bent crystal parameters
if (oe1%f_mosaic.EQ.2) then
   i = 0
   call crystal_bent_moments(oe1,i)
   call crystal_bent_ml_calc_c(oe1,psi_h,i)
endif
if (oe1%f_mosaic.EQ.3) then
   i = 0
   call crystal_bent_moments(oe1,i)
   call crystal_bent_pp_calc_beta(oe1,psi_h,i)
endif

END SUBROUTINE CRYSTAL_SET

!-------------------------------------------------------------------------------

!C+++
!C       SUBROUTINE        CRYSTAL_SET_VECTORS
!C
!C       PURPOSE           Calculates useful vectors of the crystal
!C                         (Bragg direction corrected for refraction)
!C
!C       MODIFIED          M. Sanchez del Rio
!C 
!C                         2012-12-07  created
!C
!C---

SUBROUTINE CRYSTAL_SET_VECTORS (oe1)

implicit none
                               
type(oeSetup),intent(inout)   :: oe1
real(kind=skr),dimension(3)   :: VNOR,VPAR,AXIS,BH ,BH0,vtmp
real(kind=skr),dimension(3)   :: VIN_BRAGG,VIN_BRAGG_UNCORR,VIN_BRAGG_ENERGY
real(kind=skr),dimension(3)   :: VOUT_BRAGG,VOUT_BRAGG_UNCORR,VOUT_BRAGG_ENERGY

integer(kind=ski) :: i_debug=0

!------------------------------------------------------------------------------
!C initialize the vectors  (vvin at the exact Bragg angle)
!C

! normal to the surface always on Z axis
        vnor(1) = 0.0D0
        vnor(2) = 0.0D0
        vnor(3) = 1.0D0
! vector perpendicular to vnor for defining the diffraction plane
! (diffraction plane contains vnor and vpar)
        vpar(1) = 0.0D0
        vpar(2) = 1.0D0
        vpar(3) = 0.0D0
 
! compute a vector perpendicular to both vnor and vpar
!bug fixed srio 20130206: follow right hand ref
!call cross(vnor,vpar,axis)
        call cross(vpar, vnor, axis)
        call norm(axis,axis)

! normal to Bragg planes obtaining by rotating vnor an angle equal to 
! minuns asymmetry angle (-alphaXOP) around X using 
! rodrigues rotation (in the screw direction (cw) when looking in 
! the axis direction)
!------------  
        call RODRIGUES(vnor, axis,-oe1%a_bragg,bh0)
        if (i_debug .gt. 0) then
          write(i_debug,*) 'CRYSTAL_SET: alpha_xop: ',  oe1%a_bragg*180/pi
          write(i_debug,*) 'CRYSTAL_SET: alpha: ',180.0-oe1%a_bragg*180/pi
          write(i_debug,*) 'CRYSTAL_SET: chi: ',   90.0-oe1%a_bragg*180/pi
          write(i_debug,*) 'CRYSTAL_SET: vnor: ',vnor
          write(i_debug,*) 'CRYSTAL_SET: axis: ',axis
          write(i_debug,*) 'CRYSTAL_SET: bh0: ', bh0
        endif
        bh = bh0 * TWOPI/oe1%crystalData%d_spacing

222 continue

!
! Bragg position (uncorrected)
!
        call rodrigues(-bh0, axis, (pi/2-abs(oe1%graze)) , vin_bragg_uncorr)  


        if ((dot_product(bh0,vin_bragg_uncorr)) .gt. 0) then 
           print*,'CRYSTAL_SET: H.vin: ',dot_product(bh0,vin_bragg_uncorr)
           print*,'CRYSTAL_SET: H.k_in > 0   => aborted '
           stop
        endif

!
! Bragg position (corrected)
!
        call rodrigues(-bh0,axis,(pi/2-abs(oe1%theta_b_cor)),vin_bragg) 

!
! Bragg position (for energy scan)
!
        call rodrigues(-bh0,axis,(pi/2-abs(oe1%theta_b_cor)),vin_bragg_energy) 


!
! be sure that we always enter from the upper side 
! Bug reported by sonia.francoual at desy.de because it produced
! crystal reflectivities larger than one.
!
        if (vin_bragg_uncorr(3) .gt. 0) then
            call rodrigues(-bh0, axis, -1.0*(pi/2-abs(oe1%graze)) , vin_bragg_uncorr)  
            call rodrigues(-bh0, axis, -1.0*(pi/2-abs(oe1%theta_b_cor)),vin_bragg) 
            call rodrigues(-bh0, axis, -1.0*(pi/2-abs(oe1%theta_b_cor)),vin_bragg_energy) 
        endif

!
! output directions
!
        !ewald equation for uncorrected position
        vout_bragg_uncorr = bh + vin_bragg_uncorr*oe1%q_phot_set
        vout_bragg_uncorr = vout_bragg_uncorr/ &
                    sqrt(dot_product(vout_bragg_uncorr,vout_bragg_uncorr))

        !scat direction for the others
        call scat_direction (vin_bragg,bh,vnor,oe1%q_phot_set,vout_bragg)
        vout_bragg = vout_bragg/sqrt(dot_product(vout_bragg,vout_bragg))

        call scat_direction (vin_bragg_energy,bh,vnor, &
                    oe1%q_phot_set,vout_bragg_energy)
        vout_bragg_energy = vout_bragg_energy/ &
                    sqrt(dot_product(vout_bragg_energy,vout_bragg_energy))
!
! store data
!
        oe1%VNOR = VNOR
        oe1%VPAR = VPAR
        oe1%AXIS = AXIS
        oe1%BH   = BH
        oe1%BH0  = BH0
        oe1%VIN_BRAGG         = VIN_BRAGG
        oe1%VIN_BRAGG_UNCORR  = VIN_BRAGG_UNCORR
        oe1%VIN_BRAGG_ENERGY  = VIN_BRAGG_ENERGY
        oe1%VOUT_BRAGG        = VOUT_BRAGG
        oe1%VOUT_BRAGG_UNCORR = VOUT_BRAGG_UNCORR
        oe1%VOUT_BRAGG_ENERGY = VOUT_BRAGG_ENERGY

END SUBROUTINE CRYSTAL_SET_VECTORS

!-------------------------------------------------------------------------------
! this is not used. It was a test for checking compatibility with 
! our crystal_perfect()
!-------------------------------------------------------------------------------
SUBROUTINE CRYSTAL_PERFECT_VETTIER ( &
       Q_PHOT, VIN, VOUT, &      ! inputs (ray)
       BH, SURFNOR, oe1, &       ! inputs (crystal)
       PSI_0,PSI_H,PSI_HBAR, &   ! inputs (struct fact)
       R_S, R_P,PHASE_S, PHASE_P)! output (reflectivity)

implicit none

real(kind=skr),intent(in)   ::  Q_PHOT
real(kind=skr),intent(in)   ::  VIN(3), VOUT(3), BH(3), SURFNOR(3)
type(oeSetup),intent(in)    ::  oe1
complex(kind=skx),intent(in):: psi_0,psi_h,psi_hbar

real(kind=skr),intent(out)  :: R_S, R_P,PHASE_S, PHASE_P
real(kind=skr)  :: vtemp(3)

!C****************************************************************
!C* calculates complex amplitudes (reflected or transmitted)     *
!C* from a perfect crystal in either the Bragg or Laue geometry. *
!C* uses the first expression in Zachariassen 3.130 and 3.137    *
!C****************************************************************
!C
        REAL*8 THET,LAMBDA
        REAL*8 SINTB,AKPA,RRP,RIP
        REAL*8 THETAB,B,K,GAMMA0,GAMMAH,T0

        COMPLEX*16 ZS,ZP
        COMPLEX*16 Z,ZQ,Z1,CARG1,CARG2,CDEL1,CDEL2
        COMPLEX*16 CP1,CP2,CX1,CX2,C1,C2,CDEN
        COMPLEX*16 CI
        COMPLEX*16 PSI0,PSIH,PSIHB

real(kind=skr)   :: sin_brg,sin_brg0
real(kind=skr)   :: gamma_0,gamma_h,r_lam0,sin_gra !,phot,graze
integer(kind=ski)   :: ig,i_debug=0

!---------------
!match variables
!---------------
CI = (0.0d0,1.0d0) 
psi0 = psi_0
psih = psi_h
psihb = psi_hbar

R_LAM0 = TWOPI/Q_PHOT    ! wavelength in cm
SIN_GRA = R_LAM0/oe1%crystalData%d_spacing/2.0D0  
!GRAZE = ASIN(SIN_GRA)            ! uncorrected Bragg angle

call dot (vin,surfnor,gamma_0)     ! gamma0
call dot (vout,surfnor,gamma_h)    ! gammah

call norm (bh,vtemp)
call dot (vin,vtemp,sin_brg)                   ! angle vin with BH 
call dot (oe1%vin_bragg_uncorr,vtemp,sin_brg0) ! angle vin_bragg_uncorr with BH 

!
!   zac_alpha =-((R_LAM0/oe1%crystalData%D_SPACING)**2+ &
!           2*R_LAM0* sin_brg/oe1%crystalData%D_SPACING)
!
!Call dot (bh,bh,mod2_bh)
!Call dot (vin,bh,tmp)
!zac_alpha =   (mod2_bh/q_phot**2) +  (2/q_phot) *tmp


sintb = abs(sin_gra)
!watch out this change of sign!!
gamma0 = -gamma_0
gammaH = -gamma_H

b = gamma0/gammaH
T0 = oe1%thickness 
LAMBDA  = TWOPI/Q_PHOT    ! wavelength in cm
thet = abs(asin(sin_brg))
Z=(1.D0-B)*PSI0+4.D0*B*SINTB*(SINTB-DSIN(THET))
Z=Z/2.D0                ! z de Zachariassen
if (oe1%f_refrac .eq. 0) ig=1
if (oe1%f_refrac .eq. 1) ig=2
if (oe1%f_refrac .eq. 2) ig=-1
if (oe1%f_refrac .eq. 3) ig=-2
k = abs(cos(2*asin(sintb))) !pol factor

if (i_debug.gt.0) then
  write(11,*) 'sintb: ',sintb
  write(11,*) 'gamma0: ',gamma0
  write(11,*) 'gammaH: ',gammaH
  write(11,*) 'b: ',b
  write(11,*) 'T0: ',T0
  write(11,*) 'lambda: ',lambda
  write(11,*) 'sintb: ',sintb
  write(11,*) 'thetab: ',asin(sintb)*180/pi
  write(11,*) 'thet: ',thet*180.0/pi
  write(11,*) 'thet-thetab: ',thet*180.0/pi-asin(sintb)*180/pi
  write(11,*) 'Z: ',Z
  write(11,*) 'psi0: ',psi0
  write(11,*) 'psih: ',psih
  write(11,*) 'psihb: ',psihb
  write(11,*) 'ig: ',ig
  write(11,*) 'pol fact K: ',k
endif

!end match variables
!C
!C==============         ! SIGMA polarization
        ZQ=B*PSIH*PSIHB
        Z1=ZQ+(Z*Z)
        Z1=CDSQRT(Z1)
!C==============
        AKPA=PI/(LAMBDA*GAMMA0)
        CARG1=-Z+Z1
        CARG2=-Z-Z1
        CDEL1=AKPA*(PSI0+CARG1)
        CDEL2=AKPA*(PSI0+CARG2)
!C
        CP1=-CI*CDEL1*T0
        CP2=-CI*CDEL2*T0
        CX1=CARG1/PSIH
        CX2=CARG2/PSIH
!C==============         ! stops under- and over-flows
        RRP=70.D0
        IF(DABS(DREAL(CP1)).GT.RRP) THEN
                RIP=DIMAG(CP1)
                IF(DREAL(CP1).GT.RRP) C1=CDEXP(DCMPLX(RRP,RIP))
                IF(DREAL(CP1).LT.-RRP) C1=CDEXP(DCMPLX(-RRP,RIP))
        ELSE IF (DABS(DREAL(CP1)).LE.RRP) THEN
                C1=CDEXP(CP1)
        ENDIF
        IF(DABS(DREAL(CP2)).GT.RRP) THEN
                RIP=DIMAG(CP2)
                IF(DREAL(CP2).GT.RRP) C2=CDEXP(DCMPLX(RRP,RIP))
                IF(DREAL(CP2).LT.-RRP) C2=CDEXP(DCMPLX(-RRP,RIP))
        ELSE IF (DABS(DREAL(CP2)).LE.RRP) THEN
                C2=CDEXP(CP2)
        ENDIF
!C======================
        IF (IG.GT.0) THEN               ! reflected beam
                IF(IG.EQ.1) THEN        ! Bragg case
                        CDEN=C2*CX2-C1*CX1
                        ZS=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(B))
                ELSE IF (IG.EQ.2) THEN  ! Laue case
                        CDEN=CX2-CX1
                        ZS=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(B))
                ENDIF
!C======================
        ELSE IF (IG.LT.0) THEN          ! transmitted beam
                IF(IG.EQ.-1) THEN       ! Bragg case
                        CDEN=C2*CX2-C1*CX1
                        ZS=C1*C2*(CX2-CX1)/CDEN/DSQRT(DABS(B))
                ELSE IF (IG.EQ.-2) THEN ! Laue case
                        CDEN=CX2-CX1
                        ZS=(CX2*C1-CX1*C2)/CDEN/DSQRT(DABS(B))
                ENDIF
        ENDIF
!C============================================================
!C                       ! PI polarization
        ZQ=B*K*K*PSIH*PSIHB
        Z1=ZQ+(Z*Z)                     ! q+z2 dans Zacha.
        Z1=CDSQRT(Z1)
!C==============
        AKPA=PI/(LAMBDA*GAMMA0)
        CARG1=-Z+Z1
        CARG2=-Z-Z1
        CDEL1=AKPA*(PSI0+CARG1)
        CDEL2=AKPA*(PSI0+CARG2)
!C
        CP1=-CI*CDEL1*T0
        CP2=-CI*CDEL2*T0
        CX1=CARG1/(K*PSIH)
        CX2=CARG2/(K*PSIH)
!C==============         ! stops under- and over-flows
        RRP=70.D0
        IF(DABS(DREAL(CP1)).GT.RRP) THEN
                RIP=DIMAG(CP1)
                IF(DREAL(CP1).GT.RRP) C1=CDEXP(DCMPLX(RRP,RIP))
                IF(DREAL(CP1).LT.-RRP) C1=CDEXP(DCMPLX(-RRP,RIP))
        ELSE IF (DABS(DREAL(CP1)).LE.RRP) THEN
                C1=CDEXP(CP1)
        ENDIF
        IF(DABS(DREAL(CP2)).GT.RRP) THEN
                RIP=DIMAG(CP2)
                IF(DREAL(CP2).GT.RRP) C2=CDEXP(DCMPLX(RRP,RIP))
                IF(DREAL(CP2).LT.-RRP) C2=CDEXP(DCMPLX(-RRP,RIP))
        ELSE IF (DABS(DREAL(CP2)).LE.RRP) THEN
                C2=CDEXP(CP2)
        ENDIF
!C======================
        IF (IG.GT.0) THEN               ! reflected beam
                IF(IG.EQ.1) THEN        ! Bragg case
                        CDEN=C2*CX2-C1*CX1
                        ZP=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(B))
                ELSE IF (IG.EQ.2) THEN  ! Laue case
                        CDEN=CX2-CX1
                        ZP=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(B))
                ENDIF
!C======================
        ELSE IF (IG.LT.0) THEN          ! transmitted beam
                IF(IG.EQ.-1) THEN       ! Bragg case
                        CDEN=C2*CX2-C1*CX1
                        ZP=C1*C2*(CX2-CX1)/CDEN/DSQRT(DABS(B))
                ELSE IF (IG.EQ.-2) THEN ! Laue case
                        CDEN=CX2-CX1
                        ZP=(CX2*C1-CX1*C2)/CDEN/DSQRT(DABS(B))
                ENDIF
        ENDIF
!C================================================================
write(11,*) 'cdel1: ',cdel1
write(11,*) 'cdel2: ',cdel2
write(11,*) 'cp1: ',cp1
write(11,*) 'cp2: ',cp2
write(11,*) 'cx1: ',cx1
write(11,*) 'cx2: ',cx2
write(11,*) 'c1: ',c1
write(11,*) 'c2: ',c2
write(11,*) 'zs: ',zs
write(11,*) 'zp: ',zp
        R_S = abs(ZS)
        R_P = abs(ZP)
        RETURN
END SUBROUTINE CRYSTAL_PERFECT_VETTIER


!-------------------------------------------------------------------------------
!C+++
!C       SUBROUTINE      CRYSTAL_PERFECT
!C
!C       PURPOSE         Computes the reflectivity of a perfect crystal 
!C                       according to the dynamic theory of x-ray diffraction.
!C
!C       ALGORITHM        Reference 
!C                       W. H. Zachariasen
!C                       "Theory of X-Ray diffraction in crystals"
!C                       Dover, New York, 1967.
!C
!C       MODIFIED        M. Sanchez del Rio, Feb 1996. This routine was
!C                       a part of crystal.F, and now is independent.
!C                       Aug 2012: Cleaned, documented and some variables 
!C                                 renamed.
!C                       Jan 2014: Re-Cleaned
!C
!C
!C---
!C
!C INPUT PARAMETERS:
!C   Q_PHOT:  Photon wavelength
!C   VIN:     Photon incident direction
!C   VOUT:    Photon outcoming direction
!C   BH:               Bragg planes normal
!C   SURFNOR:       Crystal surface normal
!C   OE1:     Type with crystal info, in particular:
!C            THICKNESS:       crystal thickness
!C            F_REFRAC:       Mode flag: 0:BraggDiff, 1:LaueDiff, 
!C                                       2:BraggTrans, 3:LaueTrans
!C   PSI_0:       |
!C   PSI_H:       | Electric susceptibility (proportional to structure factor)
!C   PSI_HBAR:    |
!C
!C OUTPUT PARAMETERS:
!C   R_S:       Sigma reflectivity
!C   R_P:       Pi reflectivity
!C   PHASE_S:       Sigma phase shift
!C   PHASE_P:       Pi phase shift
!C
!C
SUBROUTINE CRYSTAL_PERFECT ( &
       Q_PHOT, VIN, VOUT, &      ! inputs (ray)
       BH, SURFNOR, oe1, &       ! inputs (crystal)
       PSI_0,PSI_H,PSI_HBAR, &   ! inputs (struct fact)
       R_S, R_P,PHASE_S, PHASE_P)! output (reflectivity)

implicit none

real(kind=skr),intent(in)   ::  Q_PHOT
real(kind=skr),intent(in)   ::  VIN(3), VOUT(3), BH(3), SURFNOR(3)
type(oeSetup),intent(in)    ::  oe1
complex(kind=skx),intent(in):: psi_0,psi_h,psi_hbar

real(kind=skr),intent(out)  :: R_S, R_P,PHASE_S, PHASE_P


real(kind=skr)   :: OEXP=100.0D0   !exponential overflow limit
real(kind=skr)   :: r_lam0,phot
real(kind=skr)   :: graze,sin_gra
real(kind=skr)   :: sin_brg,sin_brg0,zac_alpha,mod2_bh,tmp
real(kind=skr)   :: gamma_0,gamma_h
real(kind=skr)   :: cry_b !,cry_t,cry_a 
real(kind=skr)   :: c_ppol,pp,qq
complex(kind=skx):: ci
complex(kind=skx):: zac_q,zac_z,ctemp,zac_x1,zac_x2
complex(kind=skx):: zac_delta1,zac_delta2,zac_c1,zac_c2
complex(kind=skx):: zac_phi1,zac_phi2
complex(kind=skx):: rcs,rcp

real(kind=skr),dimension(3)   :: vtemp,vtemp2
!
integer(kind=ski):: i_debug=0  ! set to file unit e.g., 11 (debug output)


if (i_debug.GE.1) then
  write (i_debug,*) '<><>'
  write (i_debug,*) '<><>CRYSTAL_PERFECT: ******** crystal_perfect called ********'
endif

CI        = (0.0D0,1.0D0)

PHOT    = Q_PHOT/TWOPI*TOCM  ! photon energy in eV
R_LAM0  = TWOPI/Q_PHOT       ! wavelength in cm
SIN_GRA = R_LAM0/oe1%crystalData%d_spacing/2.0D0  
GRAZE   = ASIN(SIN_GRA)      ! uncorrected Bragg angle

if (i_debug.GE.1) write(i_debug,*) '<><>CRYSTAL_PERFECT: working at energy: ',phot

!
! calculate director cosines (see Zachariasen pag 118, after formula 3.115)
!
call dot (vin,surfnor,gamma_0)     ! gamma0
call dot (vout,surfnor,gamma_h)    ! gammah
! Our crystal normal is pointing outside me medium. Zachariasen normal is
! pointing into the crystal medium (pag 112). Therefore, we must change the
! sign.
gamma_0 = -gamma_0
gamma_h = -gamma_h

if (i_debug.GE.1) then
  write(i_debug,*) '<><>CRYSTAL_PERFECT: PSI_H = ', PSI_H
  write(i_debug,*) '<><>CRYSTAL_PERFECT: PSI_HBAR = ', PSI_HBAR
  write(i_debug,*) '<><>CRYSTAL_PERFECT: PSI_0 = ', PSI_0
endif

!C
!C >>>>>>>>>>>>>>>>>>>> Perfect crystal calculation <<<<<<<<<<<<<<<<<<<
!C
!C Main calculation (symmetrical case and asym incident case)
!C I change to reflectivity formulae of Zachariasen,
!C for a definition of ETA taking into account the angle with
!C Bragg planes. MSR 6/28/90
!C

! 
! calculate alpha variable of Zachariasen (Eq 3.114b)
!
if (i_debug.GE.1) then 
   call norm (bh,vtemp)
   call dot (vin,vtemp,sin_brg)                   ! angle vin with BH 
   call dot (oe1%vin_bragg_uncorr,vtemp,sin_brg0) ! angle vin_bragg_uncorr with BH 

   !zac_alpha =-((R_LAM0/oe1%crystalData%D_SPACING)**2+ &
   !        2*R_LAM0* sin_brg/oe1%crystalData%D_SPACING)
   !write (i_debug,*) '<><>CRYSTAL_PERFECT: !!!!! zac_alpha_old: ',zac_alpha
   write (i_debug,*) '<><>CRYSTAL_PERFECT: '
   write (i_debug,*) '<><>CRYSTAL_PERFECT: theta:  ',asin(sin_brg)/PI*180
   write (i_debug,*) '<><>CRYSTAL_PERFECT: thetaB: ',asin(sin_brg0)/PI*180
   write (i_debug,*) '<><>CRYSTAL_PERFECT: graze=|thetaB|: ',graze*180/PI
   write (i_debug,*) '<><>CRYSTAL_PERFECT: bh: ',bh
end if

Call dot (bh,bh,mod2_bh)
Call dot (vin,bh,tmp)
zac_alpha =   (mod2_bh/q_phot**2) +  (2/q_phot) *tmp
if (i_debug.GE.1) then 
   write (i_debug,*) '<><>CRYSTAL_PERFECT: !!!!! R_LAM0: ',R_LAM0
   write (i_debug,*) '<><>CRYSTAL_PERFECT: !!!!! D_SPACING: ',oe1%crystalData%D_SPACING
end if


!*
!* Transmission (Laue) case and general Reflection (Bragg) case
!* NB: The general Bragg case is commented (cc). It is not yet being
!* used. We still use the Thick Crystal Approximation case.
!*
!C
!C   This part has been written by G.J. Chen and M. Sanchez del Rio. 
!C   We use the formula [3.130] of Zachariasen's book.
!C


!C
!C changed b factor to its vectorial value (Zachariasen, [3.115])
!C
if (i_debug.GE.1) then
  cry_b = gamma_0/gamma_h
  write (i_debug,*) 'CRYSTAL_PERFECT: b(approx)= ',cry_b
end if
!C numerator
call dot (surfnor,vin,cry_b)
cry_b = cry_b*q_phot
!C denominator
call scalar(vin,q_phot,vtemp)
call vsum (vtemp,bh,vtemp2)
call dot (surfnor,vtemp2,tmp)
!C ratio
cry_b = cry_b / tmp
if (i_debug.GE.1) then
  write (i_debug,*) 'CRYSTAL_PERFECT: b(exact)= ',cry_b
end if

if (i_debug.GE.1) then
  write(i_debug,*) '<><>CRYSTAL_PERFECT: zac_alpha: ',zac_alpha
  ! this is the approximated alpha, eq. 3.116, with signs!!
  tmp = 2*(asin(sin_brg0)-asin(sin_brg)) *  sin(2*asin(sin_brg0))
  write(i_debug,*) '<><>CRYSTAL_PERFECT: zac_alpha approx: ', tmp
  write(i_debug,*) '<><>CRYSTAL_PERFECT: thetaB-theta=',asin(sin_brg0)-asin(sin_brg)
  Call dot (vin,bh,tmp)
  write(i_debug,*) '<><>CRYSTAL_PERFECT: Vin.BH=',tmp
endif
!C

! Zachariasen Eqs. 3.123
zac_q = cry_b*psi_h*psi_hbar      
zac_z = (1.0D0-cry_b)*0.5D0*psi_0 + cry_b*0.5D0*zac_alpha

if (i_debug.GE.1) then
  ! y of Zachariasen (not needed)
  ! this is y of Zachariasen, eq. 3.141
  !TODO: check the formula, replace  sqrt(psi_h*psi_hbar) by cdabs(psi_h) ??
  tmp = real((1.0D0-cry_b)*psi_0+cry_b*zac_alpha)/(2.0D0*sqrt(abs(cry_b)) &
        *1.0d0*sqrt(psi_h*psi_hbar))
  write(i_debug,*) '<><>CRYSTAL_PERFECT: y of Zachariasen: ',tmp
endif
!

!C
!C s-polarization. Zachariasen Eqs 3.122, 3.121, 3.126 and 3.132
!C
ctemp      = sqrt(zac_q  + zac_z**2)
zac_x1     = (-1.0d0*zac_z+ctemp)/psi_hbar
zac_x2     = (-1.0d0*zac_z-ctemp)/psi_hbar
zac_delta1 = 0.5d0*(psi_0-zac_z+ctemp)
zac_delta2 = 0.5d0*(psi_0-zac_z-ctemp)

zac_phi1 = twopi/(gamma_0)/r_lam0*zac_delta1
zac_phi2 = twopi/(gamma_0)/r_lam0*zac_delta2
zac_c1   = -1.d0*ci*oe1%thickness*zac_phi1
zac_c2   = -1.d0*ci*oe1%thickness*zac_phi2

if (i_debug.GE.1) then
  write(i_debug,*) '<><>CRYSTAL_PERFECT: ctemp: ',ctemp
endif
!C
!C a very big exponential produces numerical overflow. If so, the value
!C is changed artificially to avoid the overflow. This is equivalent to 
!C use the thick crystal approximation
!C
!C       oexp = 100.d0
if (Dreal(zac_c1).gt.oexp.or.Dreal(zac_c2).gt.oexp) then
 if (i_debug.GE.1) write(*,*) 'CRYSTAL: exponential overflow. Corrected.'
 write(i_debug,*) 'CRYSTAL s-pol: exponential overflow. Corrected.'
 if (Dreal(zac_c1).gt.oexp) zac_c1 = oexp+ci*Dimag(zac_c1)
 if (Dreal(zac_c1).lt.-oexp) zac_c1 = -oexp+ci*Dimag(zac_c1)
 if (Dreal(zac_c2).gt.oexp) zac_c2 = oexp+ci*Dimag(zac_c2)
 if (Dreal(zac_c2).lt.-oexp) zac_c2 = -oexp+ci*Dimag(zac_c2)
endif
!C
zac_c1 = exp(zac_c1)
zac_c2 = exp(zac_c2)


!C
if (oe1%f_refrac.eq.0) then
  rcs = zac_x1*zac_x2*(zac_c1-zac_c2)/(zac_c2*zac_x2-zac_c1*zac_x1) ! bragg D
else if (oe1%f_refrac.eq.1) then 
  rcs = zac_x1*zac_x2*(zac_c1-zac_c2)/(zac_x2-zac_x1)               ! laue D
else if (oe1%f_refrac.eq.2) then 
  rcs = zac_c1*zac_c2*(zac_x1-zac_x2)/(zac_c2*zac_x2-zac_c1*zac_x1) ! bragg T
else if (oe1%f_refrac.eq.3) then 
  rcs = (zac_x2*zac_c1-zac_x1*zac_c2)/(zac_x2-zac_x1)               ! laue T
endif
if (i_debug.GE.1) then
  write(i_debug,*) '<><>CRYSTAL_PERFECT: zac_c1: ',zac_c1
  write(i_debug,*) '<><>CRYSTAL_PERFECT: zac_c2: ',zac_c2
  write(i_debug,*) '<><>CRYSTAL_PERFECT: zac_x1: ',zac_x1
  write(i_debug,*) '<><>CRYSTAL_PERFECT: zac_x2: ',zac_x2
endif

!C
!C p-polarization
!C
! todo: check if change graze by incident angle? 
c_ppol = abs(cos(2.0d0*graze))

ctemp      = sqrt(zac_q*c_ppol**2  + zac_z**2)
zac_x1     = (-1.0d0*zac_z+ctemp)/(psi_hbar*c_ppol)
zac_x2     = (-1.0d0*zac_z-ctemp)/(psi_hbar*c_ppol)
zac_delta1 = 0.5d0*(psi_0-zac_z+ctemp)
zac_delta2 = 0.5d0*(psi_0-zac_z-ctemp)
zac_phi1   = twopi/(gamma_0)/r_lam0*zac_delta1
zac_phi2   = twopi/(gamma_0)/r_lam0*zac_delta2
zac_c1     = -1.0d0*ci*oe1%thickness*zac_phi1
zac_c2     = -1.0d0*ci*oe1%thickness*zac_phi2
!C
!C a very big exponential produces numerical overflow. If so the value
!C is changed to avoid the overflow. This is equivalent to the thick
!C crystal approximation
!C
if (Dreal(zac_c1).gt.oexp.or.Dreal(zac_c2).gt.oexp) then
  if (i_debug.GE.1) write(*,*) 'CRYSTAL: exponential overflow. Corrected.'
  write(i_debug,*) 'CRYSTAL p-pol: exponential overflow. Corrected.'
  if (Dreal(zac_c1).gt.oexp) zac_c1 = oexp+ci*Dimag(zac_c1)
  if (Dreal(zac_c1).lt.-oexp) zac_c1 = -oexp+ci*Dimag(zac_c1)
  if (Dreal(zac_c2).gt.oexp) zac_c2 = oexp+ci*Dimag(zac_c2)
  if (Dreal(zac_c2).lt.oexp) zac_c2 = -oexp+ci*Dimag(zac_c2)
endif
!C
zac_c1 = exp(zac_c1)
zac_c2 = exp(zac_c2)

if (oe1%f_refrac.eq.0) then
  rcp = zac_x1*zac_x2*(zac_c1-zac_c2)/(zac_c2*zac_x2-zac_c1*zac_x1) ! bragg D
else if (oe1%f_refrac.eq.1) then 
  rcp = zac_x1*zac_x2*(zac_c1-zac_c2)/(zac_x2-zac_x1)             ! laue D
else if (oe1%f_refrac.eq.2) then 
  rcp = zac_c1*zac_c2*(zac_x1-zac_x2)/(zac_c2*zac_x2-zac_c1*zac_x1) ! bragg T
else if (oe1%f_refrac.eq.3) then 
  rcp = (zac_x2*zac_c1-zac_x1*zac_c2)/(zac_x2-zac_x1)             ! laue T
endif

! note division by |b| in intensity (thus sqrt(|b|) in amplitude) 
! for power balance (see Zachariasen pag. 122)
!
! this factor only applies to diffracted beam, not to transmitted beams
! changed srio@esrf.eu 20130131, see private communication J. Sutter (DLS)
if (oe1%f_refrac.le.1) then
  rcs = (1.0d0/sqrt(abs(cry_b)))*rcs
  rcp = (1.0d0/sqrt(abs(cry_b)))*rcp
endif

if (i_debug.GE.1) then
  write(i_debug,*) '<><>CRYSTAL_PERFECT: rcs: ',rcs
  write(i_debug,*) '<><>CRYSTAL_PERFECT: rcp: ',rcp
  write(i_debug*10,*) tmp,zac_alpha, &
        (asin(sin_brg)-asin(sin_brg0))*180/pi,&
        (asin(abs(sin_brg))-asin(abs(sin_brg0)))*180/pi,(abs(rcs))**2
endif

!
! compute moduli and phases
!
R_S  = ABS(RCS)
PP   = DREAL(RCS)
QQ   = DIMAG(RCS)
CALL ATAN_2        (QQ,PP,PHASE_S)
R_P  = ABS(RCP)
PP   = DREAL(RCP)
QQ   = DIMAG(RCP)
CALL ATAN_2        (QQ,PP,PHASE_P)

RETURN
END SUBROUTINE CRYSTAL_PERFECT



!C+++
!C       SUBROUTINE        CRYSTAL_MOSAIC
!C
!C       PURPOSE           Computes the reflectivity of a mosaic crystal 
!C                         according to the dynamic theory of x-ray diffraction.
!C
!C       ALGORITHM       Reference 
!C                       W. H. Zachariasen
!C                       "Theory of X-Ray diffraction in crystals"
!C                       Dover, New York, 1967.
!C                       and
!C                       Bacon and Lowde, Acta Crystall. 1 pag 303 (1948) for 17
!C
!C       MODIFIED        M. Sanchez del Rio, Feb 1996. This routine was
!C                       a part of crystal.F, and now is independent.
!C
!C---
!C
!C INPUT PARAMETERS:
!C   Q_PHOT: Photon wavelength = 2 pi / lambda [cm^-1]
!C   VIN:     Photon incident direction
!C   BH:               Bragg planes normal *  2 pi/d [cm^-1]
!C   SURFNOR:       Crystal surface normal
!C   OE1:           type with crystal data, in particular: 
!C                  THICKNESS: crystal thickness
!C                  F_REFRAC:  Mode flag: 0:BraggDiff, 1:LaueDiff, 
!C                                        2:BraggTrans, 3:LaueTrans
!C                  SPREAD_MOS: mosaicity
!C   Q_MOS:        Q variable = pi**2*abs(psi_h*psi_hbar)/r_lam0/sin(2*graze)
!C   REFRAC:       Refraction index
!C
!C OUTPUT PARAMETERS:
!C   R_S:       Sigma reflectivity
!C   R_P:       Pi reflectivity
!C   PHASE_S:       Sigma phase shift
!C   PHASE_P:       Pi phase shift
!C   L_EXT_S:       Scondary extinction length (s-pol)
!C   L_EXT_P:       Scondary extinction length (p-pol)
!C
!C

SUBROUTINE CRYSTAL_MOSAIC ( Q_PHOT, VIN, &       ! inputs (ray)
                            BH, SURFNOR, oe1, &  ! inputs (crystal)
                            Q_MOS,REFRAC, &      ! inputs (struct fact)
   R_S, R_P,PHASE_S, PHASE_P, L_EXT_S, L_EXT_P)  ! output (reflectivity)

implicit none
                               
real(kind=skr),intent(in)   ::   Q_PHOT     !Arguments
type(oeSetup),intent(in)    ::   oe1
real(kind=skr),dimension(3),intent(in) :: VIN, BH, SURFNOR
real(kind=skr),intent(in)   ::    q_mos
Complex(kind=skx),intent(in)::    refrac

real(kind=skr),intent(out)   ::   R_S, R_P,PHASE_S, PHASE_P
real(kind=skr),intent(out)   ::   L_EXT_S, L_EXT_P

real(kind=skr),dimension(3)  ::   vtemp(3)
real(kind=skr)   ::  phot,r_lam0
real(kind=skr)   ::  graze,sin_gra
real(kind=skr)   ::  absorp    
real(kind=skr)   ::  sin_brg
real(kind=skr)   ::  sin_q_ang !,sin_q_ref
real(kind=skr)   ::  sigma_gamma,biga_mosaic,ep,omega,kpolsqrt
real(kind=skr)   ::  smallas_mosaic,smallap_mosaic,rs_mosaic
real(kind=skr)   ::  rp_mosaic


!C
!C >>>>>>>>>>>>>>>>>>>> Mosaic crystal calculation <<<<<<<<<<<<<<<<<<<<<<
!C
PHOT    = Q_PHOT/TWOPI*TOCM
R_LAM0  = TWOPI/Q_PHOT
SIN_GRA = R_LAM0/oe1%crystalData%d_spacing/2.0D0
GRAZE = ASIN(SIN_GRA)

ABSORP  = 2.0d0 * TWOPI *(-DIMAG(REFRAC)) / R_LAM0

Call norm (bh,vtemp)
Call dot (vin,vtemp,sin_brg)        
Call dot (vin,surfnor,sin_q_ang)
!C       call dot (vout,surfnor,sin_q_ref)

EP      = abs(ASIN(SIN_BRG)) - graze ! scalar def

BIGA_MOSAIC  = -1.0d0*oe1%THICKNESS*ABSORP/SIN_Q_ANG


OMEGA =(DEXP(-EP**2/2.0D0/(oe1%SPREAD_MOS)**2)) /SQRT(TWOPI)/oe1%SPREAD_MOS
!CCC       Q_MOS = pi**2*cdabs(psi_h*psi_hbar)/r_lam0/sin(2*graze)
SIGMA_GAMMA =  OMEGA*Q_MOS
KPOLSQRT = (cos(2*graze))**2

SMALLAS_MOSAIC = SIGMA_GAMMA/absorp
SMALLAP_MOSAIC = SMALLAS_MOSAIC*KPOLSQRT

!*
!* Transmission case
!*
if (oe1%f_refrac.eq.1) then
 rs_mosaic = sinh(smallas_mosaic*biga_mosaic) * exp(-biga_mosaic*(1+smallas_mosaic))
 rp_mosaic = sinh(smallap_mosaic*biga_mosaic) * exp(-biga_mosaic*(1+smallap_mosaic))
!*
!* Reflection case
!*
else  if (oe1%f_refrac.eq.0) then
  RS_MOSAIC = 1+SMALLAS_MOSAIC+(SQRT(1+2*SMALLAS_MOSAIC))/ DTANH(bigA_MOSAIC*SQRT(1+2*SMALLAS_MOSAIC))
  RP_MOSAIC = 1+SMALLAP_MOSAIC+(SQRT(1+2*SMALLAP_MOSAIC))/ DTANH(bigA_MOSAIC*SQRT(1+2*SMALLAP_MOSAIC))
  RS_MOSAIC = SMALLAS_MOSAIC / RS_MOSAIC
  RP_MOSAIC = SMALLAP_MOSAIC / RP_MOSAIC
!*
!* other (not implemented) cases
!*
else
 write(*,*) 'CRYSTAL_MOSAIC: Error: Mode (f_refrac) not implemented.'
 stop
endif
R_S     = SQRT(RS_MOSAIC)
R_P     = SQRT(RP_MOSAIC)
!*
!* Mean value of depht into the crystal. To be used in MIRROR
!*
L_EXT_S = 1.0D0 /SIGMA_GAMMA
L_EXT_P = 1.0D0 /SIGMA_GAMMA/KPOLSQRT
!*
!*No phase change are introduced by now. (The formulae of reflectivity 
!*are already intensity, and no complex coefficient are considered).
!*This is not important because a mosaic crystal breaks always coherence
!*
PHASE_S = 0.0D0
PHASE_P = 0.0D0
RETURN
END SUBROUTINE CRYSTAL_MOSAIC



!C+++
!C       SUBROUTINE      CRYSTAL_BENT_ML
!C
!C       PURPOSE         Computes the reflectivity of a bent crystal 
!C                       according to the dynamic theory of x-ray diffraction.
!C
!C       ALGORITHM        Reference:
!C                       
!C                        M Sanchez del Rio, C Ferrero, V Mocella (1997)  
!C                        Computer simulation of bent perfect crystal 
!C                        diffraction profiles   
!C                        Proceedings of the SPIE vol.3151: 312-323  
!C                        (and references therein)
!C
!C       MODIFIED        Nicolas Perez Bocanegra 2012
!C                       Manuel Sanchez del Rio 2014
!C
!C
!C---
!C
!C INPUT PARAMETERS:
!C   Q_PHOT:  Photon wavelength
!C   VIN:     Photon incident direction
!C   VOUT:    Photon outcoming direction
!C   BH:               Bragg planes normal
!C   SURFNOR:       Crystal surface normal
!C   OE1:     Type with crystal info, in particular:
!C            THICKNESS:       crystal thickness
!C            F_REFRAC:       Mode flag: 0:BraggDiff, 1:LaueDiff, 
!C                                       2:BraggTrans, 3:LaueTrans
!C   PSI_0:       |
!C   PSI_H:       | Electric susceptibility (proportional to structure factor)
!C   PSI_HBAR:    |
!C	REFRACT
!C
!C OUTPUT PARAMETERS:
!C   R_S:       Sigma reflectivity
!C   R_P:       Pi reflectivity
!C   PHASE_S:       Sigma phase shift (zero, not calculated)
!C   PHASE_P:       Pi phase shift (zero, not calculated)
!C
!C



!
SUBROUTINE CRYSTAL_BENT_ML (Q_PHOT, VIN, VOUT,  &  ! inputs (ray)
       BH, SURFNOR, oe1, &                         ! inputs (crystal)
       PSI_0,PSI_H,PSI_HBAR,REFRAC, STRUCT, FH,FH_BAR, & ! inputs (struct fact)
       R_S, R_P,PHASE_S, PHASE_P)                        ! output (reflectivity)

implicit none

real(kind=skr),intent(in)   ::  Q_PHOT
real(kind=skr),intent(in)   ::  VIN(3), VOUT(3), BH(3), SURFNOR(3)
type(oeSetup),intent(in)    ::  oe1
complex(kind=skx),intent(in):: psi_0,psi_h,psi_hbar
complex(kind=skx),intent(in):: refrac,fh,fh_bar,struct

real(kind=skr),intent(out)  :: R_S, R_P,PHASE_S, PHASE_P


real(kind=skr)   :: r_lam0,phot 
real(kind=skr)   :: graze 
real(kind=skr)   :: sin_brg,zac_alpha,mod2_bh,tmp
real(kind=skr)   :: gamma_0,gamma_h
real(kind=skr)   :: cry_b,cry_t,cry_a

Real(kind=skr)   :: mu, skh,sk0,deltay,y
Real(kind=skr)   :: polarFactor 

Real(kind=skr)   :: cklein

real(kind=skr),dimension(3)   :: vtemp,vtemp2
Real(kind=skr)   :: asmall,g,prodt,trod,tsmall, tsmall_old
Real(kind=skr)   :: ztmp,zr,zi,zabsq,qbig,w,v,RC1
Real(kind=skr)   :: d1,d2,d3,d4,ufr,uex,rsmall,ml_ref,fact
complex(kind=8)  :: qsmall,z,usmall,zu
integer :: i,nslab
integer(kind=ski)   :: i_polarFactor
integer(kind=ski)   :: i_debug=0

!
! retrieve basic parameters
!
PHOT   = Q_PHOT/TWOPI*TOCM               !photon energy in eV
R_LAM0 = TWOPI/Q_PHOT                    !wavelength in cm
mu     = 2.0D0*twopi*(-imag(refrac))/(r_lam0)  !absorption coeff in cm^-1

!bragg angle in rad
graze = ASIN(R_LAM0/oe1%crystalData%d_spacing/2.0D0)

! angle between incident ray and Bragg planes 
Call norm (bh,vtemp)
Call dot (vin,vtemp,sin_brg)
Call dot (surfnor,vin,gamma_0)
Call dot (surfnor,vout,gamma_h)

! debugging
if (i_debug.ge.1) then 
  write(i_debug,*) 'sqrt(psi_h*psi_hbar): & ',sqrt(psi_h*psi_hbar)
  write(i_debug,*) '\\ mu: & ',mu
  write(i_debug,*) '\\ gamma\_0 gamma\_h b\_ratio: & ',gamma_0,gamma_h,gamma_0/gamma_h
endif

! calculate alpha of Zachariasen (Eq 3.114b Zachariasen)
Call dot (bh,bh,mod2_bh)
Call dot (vin,bh,zac_alpha)
zac_alpha =   (mod2_bh/q_phot**2) +  (2/q_phot)*zac_alpha

! asymmetry b factor vectorial value (Zachariasen, [3.115])
!C numerator
call dot (surfnor,vin,cry_b)
cry_b = cry_b*q_phot
!C denominator
call scalar(vin,q_phot,vtemp)
call vsum (vtemp,bh,vtemp2)
call dot (surfnor,vtemp2,tmp)
!C ratio
cry_b = cry_b / tmp

!
! this is the y-shift between lamellae (different for Bragg and laue)
!
if (oe1%f_refrac.eq.0) then ! bragg
    deltay=2.0D0
else if (oe1%f_refrac.eq.1) then ! laue
    deltay=pi/2.0D0
endif
!
! introduce the direction of the lamella rotation. 
! if concave (R>0) the y of the successive lamellae increases
! if convex (R<0) the y of the successive lamellae reduces
deltay = sign(deltay,oe1%bent_rm)


!
! loop over polarization factor
!
do i_polarFactor=0,1
    !polarization factor
    if (i_polarfactor.eq.0) then
      polarFactor = 1.00D0  !s polarized
    else if (i_polarFactor.eq.1) then
      polarFactor = abs(cos(2.00d0*graze))!p polarized
    endif
    ! debugging
    if (i_debug.ge.1) then 
      write(i_debug,*)'Pol factor polarFactor: & ',polarFactor
    endif
    
    cklein = oe1%ml_c/(polarFactor**2)
    
    if (i_debug.gt.1) then 
       write(i_debug,*) '\\ c used for calculations: & ',cklein
    endif
    
    !TODO:  check sign, removed abs() 
    cry_a = deltay/abs(cklein)
    
    ! cry_t is the slab thickness, Delta T, same plave in Xianbo's paper
    cry_t = cry_a*((oe1%l_ext_s)/polarFactor)
    
    ! paths inside the lamella Lo and Lh (just after eq. 13)
    skh=cry_t/abs(gamma_h)
    sk0=cry_t/abs(gamma_0)
    
    
    !lamella cannot be thicker than the crystal itself
    if (cry_t.gt.oe1%thickness) cry_t = oe1%thickness 
    
    !number of slabs
    nslab=int(oe1%thickness/cry_t)
    
    !
    ! prepare variables for calculating reflectivity and transmittivity
    ! (using Eqs 7 and 8 in Caciuffo et al.)
    !
    qsmall = cry_b*psi_h*psi_hbar*polarFactor*polarFactor
    asmall = pi*cry_t/(r_lam0*abs(gamma_0))
    !g = (1.00d0-cry_b)*(-imag(psi_0))/(2.00d0*polarFactor*abs(cry_b)**0.50d0* &
    !    sqrt(abs(psi_h*psi_hbar)))
    g = (1.00d0-cry_b)*imag(psi_0)/(2.00d0*polarFactor*abs(cry_b)**0.50d0* &
        sqrt(abs(psi_h*psi_hbar)))
    
    ! initialize variables
    ml_ref = 0.00D0
    prodt = 1.00D0
    tsmall_old = 1.00D0
    
    !TODO: check the formula, replace  sqrt(psi_h*psi_hbar) by cdabs(psi_h) ??
    !y = real((1.0D0-cry_b)*(-psi_0)+cry_b*zac_alpha)/ &
    y = real((1.0D0-cry_b)*psi_0+cry_b*zac_alpha)/ &
             (2.0D0*sqrt(abs(cry_b))*polarFactor*sqrt(psi_h*psi_hbar))
    
    !kappa = imag(psi_h)/real(psi_h)
    
    if (i_debug.ge.1) then 
       write(i_debug,*) '\\ slab thickness [cm]: & ',cry_t
       write(i_debug,*) '\\ number of slabs: & ',nslab
       write(i_debug,*) '\\ y Zac: & ',y
       write(i_debug,*) '\\ alpha Zac: & ',zac_alpha
       write(i_debug,*) '\\ g: & ',g
       !write(i_debug,*) '\\ kappa: & ',kappa
       !write(i_debug,*) '\\ peak reflectivity: & ', &
       !                 1.0+2*g*g-2.0*sqrt(abs(g*g*(1+g*g)-kappa*kappa))
       write(i_debug,*) '\\ IntRefMos: & ',pi*pi*cdabs(psi_h)/(2.0*mu*r_lam0)
       write(i_debug,*) ' '
    endif
    
    ! loop over the slabs
    DO i=1,nslab
            ! calculate reflectivity of the single slab for the given y
            ! (y enters in z, etc.)
            ! TODO: use crystal_perfect() instead
            ztmp = polarFactor*sqrt(abs(psi_h*psi_hbar))*sqrt(abs(cry_b))
            zr = y*ztmp                                                     
            zi = g*ztmp
            z = cmplx(zr,zi)
            zabsq = abs(z)**2  
            qbig = abs(qsmall+z*z)
            usmall = sqrt(qsmall+z*z)
            w = imag(usmall)
            v = real(usmall)                                                   
            zu = conjg(z)*usmall
    
            ! old 
            !RC1 = exp(2.0*cry_a*(1.0+cry_b)*g/(1.0-cry_b))
            ! new (srio)
            RC1 = exp(-2.0*asmall*imag(z)*(cry_b+1.0)/(cry_b-1.0))
           
            ufr=0.0D0
            uex=0.0D0
            tsmall=0.0D0
            rsmall=0.0D0
            d1=qbig+(qbig+zabsq)*abs(sinh(asmall*w))**2
            d2=(qbig-zabsq)*dsin(asmall*v)**2
            d3=real(-zu)*dsinh(2.0D0*asmall*w)
            d4=imag(zu)*dsin(2.0D0*asmall*v)
            If(oe1%f_refrac.eq.0)then ! bragg
                ufr=(sin(asmall*v)**2+abs(dsinh(asmall*w))**2)* &
                    abs(cry_b)*abs(psi_hbar)**2
                uex=RC1
                rsmall= polarFactor*polarFactor*ufr/(d1-d2+d3+d4)
                tsmall =qbig*uex/(d1-d2+d3+d4)
                ! TODO: alternative way (VERY SLOW...)
                !call CRYSTAL_BENT_ML_ZAC( & 
                !      oe1%f_refrac,y,cry_b, r_lam0, cry_t, gamma_0, &
                !   polarFactor,psi_0,psi_h,psi_hbar,rsmall,tsmall)
                prodt=prodt*tsmall_old
                ml_ref=ml_ref+(rsmall*prodt*exp(-1.0d0*mu*skh*(i-1)))
            Else If(oe1%f_refrac.eq.1)then  ! laue
                ! laue (see formula in Xianbo's paper eq. 13
                ! (note that this formula for Laue is different
                !  than for Bragg)
                uex=RC1
                ufr = (sin(asmall*v)**2+abs(sinh(asmall*w))**2)*abs(cry_b)* &
                      abs(psi_h*psi_hbar)
                !rsmall = polarFactor*polarFactor*ufr/qbig*exp(-mu*0.5D0*(sk0+skh))
                rsmall = polarFactor*polarFactor*ufr/qbig 
                !rsmall = rsmall/abs(cry_b)
    
                tsmall = (d1-d2-d3-d4)*uex/qbig !*exp(-mu*0.5D0*(sk0+skh)) !t_i
                ! TODO: alternative way (VERY SLOW...)
                !call CRYSTAL_BENT_ML_ZAC(  &
                !      oe1%f_refrac,y,cry_b, r_lam0, cry_t, gamma_0, &
                !  polarFactor,psi_0,psi_h,psi_hbar,rsmall,tsmall)
    
                prodt=prodt*tsmall_old
                ml_ref = ml_ref+(rsmall*prodt*exp(-1.0D0*mu*skh*(nslab-i)))
           endIF
     
           tsmall_old=tsmall  !t_i-1
           y = y+deltay
    enddo
    
    !
    ! prepare outputs
    !
    if (i_polarFactor .eq. 0) then 
        R_S=sqrt(ml_ref)  ! to return amplitude
    else if (i_polarFactor .eq. 1) then 
        R_P=sqrt(ml_ref)  ! to return amplitude
    end if
        
end do ! i_polarFactor

!this method does not give phase changes
PHASE_S = 0.0
PHASE_P = 0.0

return

END SUBROUTINE CRYSTAL_BENT_ML

!
! Auxiliar routine (private) to be used by CRYSTAL_BENT_ML 
! It implements the reflectivity and transmitivity of a single crystallite
! using the Zachariasen theory
!
! NOT USED, AS IT SLOWS THE RUNNING TIME...
! TODO: remove it and use CRYSTAL_PERFECT instead
!
!

SUBROUTINE CRYSTAL_BENT_ML_ZAC( braggOrLaue,y,cry_b, r_lam0, cry_t, gamma_0, &
           polarFactor,psi_0,psi_h,psi_hbar,rsmall,tsmall)

implicit none

real(kind=skr),intent(in)  :: cry_b,cry_t,gamma_0,polarFactor,y,r_lam0
complex(kind=skx),intent(in):: psi_0,psi_h,psi_hbar
integer(kind=ski)  :: braggOrLaue
real(kind=skr),intent(out)  :: rsmall,tsmall


Real(kind=skr)    :: asmall,g,zr,zi,zabsq,qbig,w,v,RC1
Real(kind=skr)    :: d1,d2,d3,d4,ufr,uex
complex(kind=skx) :: qsmall,z,usmall,zu, ztmp
integer :: i

!print*,'    in CRYSTAL_BENT_ML_ZAC'
!
! prepare variables for calculating reflectivity and transmittivity
! (using Eqs 7 and 8 in Caciuffo et al.)
!
qsmall = cry_b*psi_h*psi_hbar*polarFactor*polarFactor
asmall = pi*cry_t/(r_lam0*abs(gamma_0))
g = (1.00d0-cry_b)*(-imag(psi_0))/(2.00d0*polarFactor*abs(cry_b)**0.50d0* &
    sqrt(abs(psi_h*psi_hbar)))

! calculate reflectivity of the single slab for the given y
! (y enters in z, etc.)
ztmp = polarFactor*sqrt(abs(psi_h*psi_hbar))*sqrt(abs(cry_b))
zr = y*ztmp                                                     
zi = g*ztmp
z = cmplx(zr,zi)
zabsq = abs(z)**2  
qbig = abs(qsmall+z*z)
usmall = sqrt(qsmall+z*z)
w = imag(usmall)
v = real(usmall)                                                   
zu = conjg(z)*usmall
! old 
! RC1 = exp(2.0*cry_a*(1.0+cry_b)*g/(1.0-cry_b))
! new (srio)
RC1 = exp(-2.0*asmall*imag(z)*(cry_b+1.0)/(cry_b-1.0))

ufr=0.0D0
uex=0.0D0
tsmall=0.0D0
rsmall=0.0D0
d1=qbig+(qbig+zabsq)*abs(sinh(asmall*w))**2
d2=(qbig-zabsq)*dsin(asmall*v)**2
d3=real(-zu)*dsinh(2.0D0*asmall*w)
d4=imag(zu)*dsin(2.0D0*asmall*v)
If(braggOrLaue.eq.0)then ! bragg
     ufr=(sin(asmall*v)**2+abs(dsinh(asmall*w))**2)*abs(cry_b)*abs(psi_hbar)**2
     uex=RC1
     rsmall= polarFactor*polarFactor*ufr/(d1-d2+d3+d4)
     tsmall =qbig*uex/(d1-d2+d3+d4)
Else If(braggOrLaue.eq.1)then ! laue 
     uex=RC1
     ufr = (sin(asmall*v)**2+abs(sinh(asmall*w))**2)*abs(cry_b)* &
           abs(psi_h*psi_hbar)
     ! srio: I believe the absorption due to the path inside the 
     ! crystallite should not be here. Removed
     ! rsmall = polarFactor*polarFactor*ufr/qbig*exp(-mu*0.5D0*(sk0+skh)) !r_i
     ! tsmall = (d1-d2-d3-d4)*uex/qbig*exp(-mu*0.5D0*(sk0+skh)) !t_i
     rsmall = polarFactor*polarFactor*ufr/qbig !*exp(-mu*0.5D0*(sk0+skh)) !r_i
     tsmall = (d1-d2-d3-d4)*uex/qbig !*exp(-mu*0.5D0*(sk0+skh)) !t_i
endIF

return

END SUBROUTINE CRYSTAL_BENT_ML_ZAC


!C+++
!C       SUBROUTINE      CRYSTAL_BENT_PP
!C
!C       PURPOSE         Computes the reflectivity of a bent crystal
!C                       according to the dynamic theory of x-ray diffraction.
!C
!C
!C       MODIFIED        Nicolas Perez Bocanegra 2012. 
!C
!C---
!C
!C INPUT PARAMETERS:
!C   Q_PHOT:  Photon wavelength
!C   VIN:     Photon incident direction
!C   VOUT:    Photon outcoming direction
!C   BH:               Bragg planes normal
!C   SURFNOR:       Crystal surface normal
!C   OE1:     Type with crystal info, in particular:
!C            THICKNESS:       crystal thickness
!C            F_REFRAC:       Mode flag: 0:BraggDiff, 1:LaueDiff,
!C                                       2:BraggTrans, 3:LaueTrans
!C   PSI_0:       |
!C   PSI_H:       | Electric susceptibility (proportional to structure factor)
!C   PSI_HBAR:    |
!C   REFRACT
!C
!C OUTPUT PARAMETERS:
!C   R_S:       Sigma reflectivity (amplitude)
!C   R_P:       Pi reflectivity (amplitude)
!C   PHASE_S:   Sigma phase shift (zero)
!C   PHASE_P:   Pi phase shift (zero)
!C
!C
SUBROUTINE CRYSTAL_BENT_PP ( &
            Q_PHOT, VIN, VOUT, &      ! inputs (ray)
            BH, SURFNOR, oe1, &       ! inputs (crystal)
            PSI_0,PSI_H,PSI_HBAR, &   ! inputs (struct fact)
            STRUCT, REFRAC, &         ! inputs (struct fact) todo: rm?
            R_S, R_P, PHASE_S, PHASE_P) ! output (reflectivity)

implicit none

real(kind=skr),               intent(in)   ::  Q_PHOT
real(kind=skr),dimension(3),  intent(in)   ::  VIN, VOUT, BH, SURFNOR
type(oeSetup),                intent(in)   ::  oe1
complex(kind=skx),            intent(in)   :: psi_0,psi_h,psi_hbar,refrac,struct
real(kind=skr),               intent(out)  :: R_S,R_P,PHASE_S,PHASE_P

real(kind=skr)                :: IR, IT, expon, eloss
real(kind=skr),dimension(3)   :: vtemp,vtemp2
Real(kind=skr)   :: d_spacing 
Real(kind=skr)   :: polarFactor,betac, beta,b22,d22
real(kind=skr)   :: r_lam0,phot,G 
real(kind=skr)   :: graze 
real(kind=skr)   :: zac_alpha,mod2_bh,tmp
real(kind=skr)   :: gamma_0,gamma_h
real(kind=skr)   :: cry_b,cry_t,cry_a 
real(kind=skr)   :: eps,pend_len
Real(kind=skr)   :: mu,y
Real(kind=skr)   :: IRp,IRm,ITp,ITm,cip,cim,cep,cem
integer(kind=ski)  :: i_polarFactor,i_debug=0

!
! retrieve basic parameters
!
d_spacing = oe1%crystalData%d_spacing    !d_spacing in cm
PHOT = Q_PHOT/TWOPI*TOCM                 !photon energy in eV
R_LAM0 = TWOPI/Q_PHOT                    !wavelength in cm
mu = 2.0D0*twopi*(-imag(refrac))/(r_lam0)  !absorption coeff in cm^-1

graze = oe1%graze

! angle between incident ray and Bragg planes 
Call dot (surfnor,vin,gamma_0)
Call dot (surfnor,vout,gamma_h)

! calculate alpha of Zachariasen (Eq 3.114b Zachariasen)
Call dot (bh,bh,mod2_bh)
Call dot (vin,bh,zac_alpha)
zac_alpha = (mod2_bh/q_phot**2) +  (2/q_phot)*zac_alpha

! asymmetry b factor to its vectorial value (Zachariasen, [3.115])
! numerator
Call dot (surfnor,vin,cry_b)
cry_b = cry_b*q_phot
! denominator
Call scalar(vin,q_phot,vtemp)
Call vsum (vtemp,bh,vtemp2)
Call dot (surfnor,vtemp2,tmp)
! ratio
cry_b = cry_b / tmp

if ( (abs(imag(psi_0))) .le. 1d-10) then 
   eps= 1.0  ! avoid division by zero
else
   eps= imag(sqrt(psi_h*psi_hbar))/imag(psi_0)
endif

G = oe1%pp_G

!
! this is a loop in polarization factor
!
do i_polarFactor=0,1
    
    !polarization factor
    if (i_polarFactor.eq.0) then
      polarFactor = 1.00D0  !s polarized
    else if (i_polarFactor.eq.1) then
      polarFactor = abs(cos(2.00d0*graze))!p polarized
    endif
    !pendellosung length in intensity
    pend_len = r_lam0*abs(gamma_0)/(polarFactor*abs(psi_h))/sqrt(abs(cry_b))
    ! this is beta as defined by Xianbo
    beta=G*r_lam0/(polarFactor*abs(psi_h)*sqrt(gamma_0*gamma_h))
    betac = pi/(2.0D0*pend_len) !critical strain gradient
    
    !
    ! this is the absorbed intensity by the created wavefields
    !
    if(abs(beta).gt.1d-6)then   
      !intensity fraction of the created wavefield
      eloss = exp(-2.0d0*pi*betac/abs(beta))
    else
      eloss = 0.0d0
    endif
    
    !TODO: check the formula, replace  sqrt(psi_h*psi_hbar) by cdabs(psi_h) ??
    !conversion from angular form (zac_alpha) to y-Zachariasen  Eq [3.141]
    !note that sqrt(psi_h*psi_hbar) is complex and we need a real!
    !y=real((1.0D0-cry_b)*real(psi_0)+cry_b*zac_alpha)/(2.0D0*sqrt(abs(cry_b))*polarFactor*sqrt(psi_h*psi_hbar))  
    y = real((1.0D0-cry_b)*psi_0+cry_b*zac_alpha)/ &
        (2.0D0*sqrt(abs(cry_b))*polarFactor*sqrt(psi_h*psi_hbar))  
    
    cep = y + sqrt(y*y+cry_b)
    cem = y - sqrt(y*y+cry_b)
    cip = y-beta*oe1%thickness + sqrt( (beta*oe1%thickness-y)**2+cry_b)
    cim = y-beta*oe1%thickness - sqrt( (beta*oe1%thickness-y)**2+cry_b)
    
    
    if (i_debug.ge.1) then 
      write(i_debug,*) ' '
      write(i_debug,*) '--------------- P: ',i_polarfactor,polarFactor
      write(i_debug,*) '\\ beta,betac: & ',beta,betac
      write(i_debug,*) '\\ G: & ',G
      write(i_debug,*) '\\ eloss, 1-eloss: & ',eloss,1-eloss
      write(i_debug,*) '\\ graze [deg]: & ',graze*180/pi
      write(i_debug,*) '\\ alpha Zac: & ',zac_alpha
      write(i_debug,*) '\\ y Zac: & ',y
      write(i_debug,*) '>>> c_pp: & ',oe1%l_ext_s*beta
      write(i_debug,*) '>>> beta: & ',beta
      write(i_debug,*) '>>> D theta[eta]: ',beta*oe1%thickness
      write(i_debug,*) '>>> D theta[rad]: ',(beta*oe1%thickness)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)
      write(i_debug,*) '>>> DE/E = ', (beta*oe1%thickness)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)/ tan(oe1%graze)
      write(i_debug,*) '>>> DE [eV] = ', (beta*oe1%thickness)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)/ tan(oe1%graze)*phot
    endif
    
    !
    ! calculate reflectivities 
    !
    if ((cep*cip).ge.0) then
      ITp = ( cry_b/(cip*cip+cry_b) ) * ( cry_b/(cep*cep+cry_b) )
      IRp =  ITp *cep*cep/cry_b
      
      if (abs(beta).gt.1d-6) then 
        b22 = (cry_b-1.0d0)*(cep-cip)/(2.0d0*beta*oe1%thickness)
        d22 = cry_b*eps*log(cep/cip)/(beta*oe1%thickness)
      else
        !TODO: check this formula taken from Veijo's code
        b22 = 0d0
        d22 = -2d0*cip/(1+cip*cip)
      endif
      expon = exp(-mu*oe1%thickness*(1.0d0+b22+d22)/abs(gamma_0) )
    else
      print*,'cryst_bent_pp: Zero reflectivity for branch +'
      IRp = 0.0 
      ITp = 0.0
      expon = 1.0
    endif
    
    tmp = expon
    
    if ((cem*cim).ge.0) then
      ITm = ( cry_b/(cim*cim+cry_b) ) * ( cry_b/(cem*cem+cry_b) )
      IRm =  ITm *cem*cem/cry_b
      
      if (abs(beta).gt.1d-6) then 
        b22 = (cry_b-1.0d0)*(cem-cim)/(2.0d0*beta*oe1%thickness)
        d22 = cry_b*eps*log(cem/cim)/(beta*oe1%thickness)
      else
        !TODO: check this formula taken from Veijo's code
        print*,'cryst_bent_pp: Undeformed crystal: beta=',beta
        b22 = 0d0
        d22 = -2d0*cim/(1+cim*cim)
      endif
      expon = exp(-mu*oe1%thickness*(1.0d0+b22+d22)/abs(gamma_0) )
    else
      print*,'cryst_bent_pp: Zero reflectivity for branch -'
      IRm = 0.0 
      ITm = 0.0
      expon = 1.0
    endif
    
    
    !todo: check if the asymmetry factor 1/|b| must be placed here
    !IR = sqrt( (IRp+IRm)*expon*(1d0-eloss) )
    !IT = sqrt( (ITp+ITm)*expon )
    IR = (IRp+IRm)*expon*(1d0-eloss)
    IT = (ITp+ITm)*expon
    
    !
    ! to avoid problems if parameters are not correctly set
    !
    if (oe1%f_refrac.eq.0) then                !bragg diffracted
       IR=0d0
       IT=0d0
    else if (oe1%f_refrac.eq.1) then           !laue diffracted
    else if (oe1%f_refrac.eq.2) then           !bragg transmitted
       IR=0d0
       IT=0d0
    else if (oe1%f_refrac.eq.3) then           !laue transmitted
    endif
    
    !
    ! these are the outputs (amplitude)
    !
    if (i_polarFactor .eq. 0) then 
      R_S = sqrt(IR)
      if (i_debug.ge.1) then 
          write(i_debug,*) 'i_polarFactor,R: ',i_polarFactor,R_S
      endif
    else if (i_polarFactor .eq. 1) then 
      R_P = sqrt(IR)
      if (i_debug.ge.1) then 
          write(i_debug,*) 'i_polarFactor,R: ',i_polarFactor,R_P
      endif
    end if
   
enddo ! i_polarFactor

!phases are not calculated by this method
PHASE_S = 0.0
PHASE_P = 0.0

RETURN
END SUBROUTINE CRYSTAL_BENT_PP


!C+++
!C       SUBROUTINE      PLOT_GLE
!C
!C       PURPOSE         Prepares a sketch of the geometry: input file for
!C                       the GLE application 
!C
!C       MODIFIED        M Sanchez del Rio
!C
!C---
!C
!C INPUT PARAMETERS:
!C   OE1:     Type with crystal info, in particular:
!C   v0i,vhi,hi: vectors
!C   iplotangles: set to 1 for drawing angles
!C
!C OUTPUT PARAMETERS:
!C   file: diff_patt.gle
!C
!C
SUBROUTINE PLOTGLE(v0i,vHi,Hi,oe1,iplotangles)

implicit none

real(kind=skr),dimension(3),intent(in)   ::  v0i,vHi,Hi
type(oeSetup),              intent(in)   ::  oe1
integer(kind=ski),          intent(in)   ::  iplotangles

real(kind=skr),dimension(3)              ::  v0 ,vH ,H 
real(kind=skr),dimension(2)              ::  myv0 ,myvH ,myH 
integer(kind=ski)  ::  axisH=2,axisV=3! xop


v0 = v0i/sqrt(sum(v0i**2))
vH = vHi/sqrt(sum(vHi**2))
 H =  Hi/sqrt(sum( Hi**2))

myv0 =  (/ v0(axisH),v0(axisV) /)
myvH =  (/ vH(axisH),vH(axisV) /)
myH  =  (/ H(axisH) ,H(axisV) /)
!
OPEN (45,FILE="diff_pat.gle",STATUS='unknown', FORM='FORMATTED',ERR=77)

write(45,*) '!'
write(45,*) '! gle code for creating crystal diagram'
write(45,*) '!'
write(45,*) '! written by ...'
write(45,*) ''

!
! variables
!
write(45,*) 'v0_h = ',myv0(1)
write(45,*) 'v0_v = ',myv0(2)
write(45,*) 'vH_h = ',myvh(1)
write(45,*) 'vH_v = ',myvh(2)
write(45,*) 'H_h = ',myH(1)
write(45,*) 'H_v = ',myH(2)
write(45,*) 'thetaB = ',oe1%graze*180/pi
write(45,*) 'alphaX = ',oe1%a_bragg*180/pi
write(45,*) 'H$ = "$x_2$"'
write(45,*) 'V$ = "$x_3$"'
write(45,*) 'hshift = 0.2'
write(45,*) 'vshift = 0.0'
write(45,*) ''
write(45,*) 'size 3 3'
write(45,*) 'set hei 0.2'
write(45,*) ''


write(45,*) '! draw axes'
write(45,*) 'set color black'
write(45,*) 'set lwidth 0.01'
write(45,*) 'amove 1 1'
write(45,*) 'rline 0 2 arrow end'
write(45,*) 'amove 1.1 2.65'
write(45,*) 'tex V$'
write(45,*) 'amove 0.75 2.65'
write(45,*) 'tex "$\vec{n}$"'
write(45,*) 'amove 1 1'
write(45,*) 'rline 2 0 arrow end'
write(45,*) 'amove 2.7 1.1'
write(45,*) 'tex H$'
write(45,*) 'set lwidth 0'
write(45,*) ''
write(45,*) ''


write(45,*) '! draw vectors'
write(45,*) 'set color red'
write(45,*) 'amove v0_h*(-1)+1  v0_v*(-1)+1'
write(45,*) 'tex "$\vec{k^0}$"'
write(45,*) 'rline v0_h v0_v arrow end '
if (iplotangles .gt.0) then
  write(45,*) 'amove v0_h*(-1)+1  v0_v*(-1)+1-0.4'
  write(45,*) 'tex "$\theta_B$"'
endif
write(45,*) ''

if (iplotangles .gt.0) then
  write(45,*) 'set color pink'
  write(45,*) 'amove 1 1 '
  write(45,*) 'rline v0_h v0_v arrow end '
  write(45,*) 'tex "$\vec{k^0}$"'
  write(45,*) ''
endif

write(45,*) 'set color red'
write(45,*) 'amove 1 1 '
if (iplotangles .gt.0) then
  write(45,*) 'arc 0.7 180-alphaX-thetaB 180-alphaX'
  write(45,*) 'set color pink'
  write(45,*) 'arc 0.7 -alphaX -alphaX+thetaB'
  write(45,*) 'arc 0.75 -alphaX-thetaB -alphaX'
  write(45,*) 'set color red'
endif
write(45,*) 'rline vH_h vH_v  arrow end'
write(45,*) 'rmove hshift vshift'
write(45,*) 'tex "$\vec{k^H}$"'
write(45,*) ''

write(45,*) 'set color blue'
write(45,*) 'amove 1 1 '
write(45,*) 'rline H_h H_v arrow end'
write(45,*) 'rmove hshift vshift'
write(45,*) 'tex "$\vec{H}$"'
write(45,*) ''

if (iplotangles .gt.0) then
  write(45,*) 'set color skyblue'
  write(45,*) 'amove 1 1 '
  write(45,*) 'rline -H_h -H_v arrow end'
  write(45,*) 'tex "$-\vec{H}$"'

  write(45,*) 'amove 1 1 '
  write(45,*) 'arc 0.5 -90-alphaX -thetaB-alphaX arrow end'
  write(45,*) 'amove 1.0 0.3'
  write(45,*) '!set hei 0.1'
  write(45,*) '!set texscale scale'
  write(45,*) 'tex "$90-\theta_B$"'
  write(45,*) '!set hei 0.2'

  write(45,*) ''
endif

write(45,*) 'set color green'
write(45,*) 'amove 1 1 '
write(45,*) 'rline -H_v*1.5 H_h*1.5 '
if (iplotangles .gt.0) then
  write(45,*) 'amove 1 1 '
  if ( (oe1%a_bragg) .ne. 0) write(45,*) 'arc 0.9 -alphaX 0 arrow start'
  write(45,*) 'amove 1.95 0.80'
  write(45,*) 'tex "$\alpha_{X}$"'
endif
write(45,*) 'amove 1 1 '
write(45,*) 'rline  H_v*1.5 -H_h*1.5 '
write(45,*) ''

if (iplotangles .gt.0) then
  write(45,*) 'set color springgreen'
  write(45,*) 'amove 1 1 '
  if ((oe1%f_refrac .eq. 0) .or. (oe1%f_refrac .eq. 2)) then  !bragg
    write(45,*) 'arc 1.05 0 180-alphaX arrow end'
    write(45,*) 'amove 1 2.1'
    write(45,*) 'tex "$\alpha$"'
  else !laue
    write(45,*) 'arc 1.05 0 90-alphaX arrow end'
    write(45,*) 'amove 1 2.1'
    write(45,*) 'tex "$\chi$"'
  endif
  write(45,*) ''
endif

write(45,*) '!draw crystal'
write(45,*) 'set color black'
write(45,*) 'amove 2 1'
write(45,*) 'box -2 -0.1'
write(45,*) ''

print *,"File written to disk: diff_pat.gle"
return
77 continue
print *,"Error writing file: diff_pat.gle"
END SUBROUTINE PLOTGLE


!
!
!


!C+++
!C       SUBROUTINE      PLOT_GLE
!C
!C       PURPOSE         Prepares a sketch of the geometry for XOP (IDL code)
!C
!C       MODIFIED        M Sanchez del Rio
!C
!C---
!C
!C INPUT PARAMETERS:
!C   v0i,vhi,hi: vectors
!C
!C OUTPUT PARAMETERS:
!C   file: diff_patt.xop
!C
!C

SUBROUTINE PLOTIDL(v0i,vHi,Hi)

implicit none

real(kind=skr),dimension(3),intent(in)   ::  v0i,vHi,Hi
real(kind=skr),dimension(3)              ::  v0 ,vH ,H 
real(kind=skr),dimension(2)              ::  myv0 ,myvH ,myH 
integer(kind=ski)  ::  axisH=2,axisV=3! xop


v0 = v0i/sqrt(sum(v0i**2))
vH = vHi/sqrt(sum(vHi**2))
 H =  Hi/sqrt(sum( Hi**2))

myv0 =  (/ v0(axisH),v0(axisV) /)
myvH =  (/ vH(axisH),vH(axisV) /)
myH  =  (/ H(axisH) ,H(axisV) /)
!
!

OPEN (45,FILE="diff_pat.xop",STATUS='unknown', FORM='FORMATTED',ERR=77)

write(45,*) ';'
write(45,*) '; idl/pro code for creating crystal diagram'
write(45,*) ';'
write(45,*) '; written by crystal3.f90 module'
write(45,*) ''
write(45,*) 'v0_h = ',myv0(1)
write(45,*) 'v0_v = ',myv0(2)
write(45,*) 'vH_h = ',myvh(1)
write(45,*) 'vH_v = ',myvh(2)
write(45,*) 'H_h = ',myH(1)
write(45,*) 'H_v = ',myH(2)
write(45,*) 'H = "x!I2"'
write(45,*) 'V = "x!I3"'
write(45,*) 'hshift = 0.2'
write(45,*) 'vshift = 0.0'
write(45,*) ''
write(45,*) '; draw axes'
write(45,*) 'plot,[-2,2],[-2,2],/nodata,xstyle=-1,ystyle=-1,xrange=[0,0],yrange=[0,0],tick=0.00001,charsize=0.0001'
write(45,*) 'arrow,0,0,0,2,/data'
write(45,*) 'arrow,0,0,2,0,/data'
write(45,*) 'xyouts,2.0,0.0, H, charsize=3'
write(45,*) 'xyouts,0.0,2.0, V, charsize=3'
write(45,*) ''
write(45,*) ''
write(45,*) '; draw vectors'
write(45,*) ';set color red'
write(45,*) ''
write(45,*) ''
write(45,*) 'arrow,-v0_h,-v0_v,0,0,/data,color=2'
write(45,*) 'xyouts, -v0_h+hshift,-v0_v+vshift, "k!Io",charsize=3'
write(45,*) 'arrow,0,0,vH_h,vH_v,/data,color=2'
write(45,*) 'xyouts, vH_h+hshift,vH_v+vshift, "k!IH",charsize=3'


write(45,*) ';set color blue'
write(45,*) ''

write(45,*) 'arrow,0,0,H_h,H_v,/data,color=4'
write(45,*) 'xyouts, H_h+hshift,H_v+vshift, "H",charsize=3'

write(45,*) ';set color green'
write(45,*) ''


write(45,*) 'arrow, 0,0, -H_v*1.5, H_h*1.5 ,/data,color=3,hsize=0.0'
write(45,*) 'arrow, 0,0,  H_v*1.5,-H_h*1.5 ,/data,color=3,hsize=0.0'

write(45,*) ';draw crystal'
write(45,*) ''

write(45,*) 'x1 = -0.8 & y1 = -0.1'
write(45,*) 'x2 =  0.8 & y2 =  0.0'
write(45,*) 'XC=[X1,X2,X2,X1,X1]'
write(45,*) 'YC=[Y1,Y1,Y2,Y2,Y1]'
write(45,*) 'Plots, XC, YC, color=5, /data'

print *,"File written to disk: diff_pat.xop"
return
77 continue
print *,"Error writing file: diff_pat.xop"
END SUBROUTINE PLOTIDL

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!C+++
!C       SUBROUTINE      CRYSTAL_BENT_ML_CALC_C
!C
!C       PURPOSE         Computes the "c" parameter of the ML model
!C
!C       MODIFIED        M Sanchez del Rio
!C
!C---
!C
!C INPUT PARAMETERS:
!    OE1:     Type with crystal info, in particular:
!    I_UNIT:  >0 write output in unit=i_unit
!    psi_h: the structure factor
!    i_unit: if >0, writes info in file unit i_unit
! 
!  OUTPUT PARAMETERS:
!    sets the following parameters
!        oe1%bent_ml_c
! 
!  ALGORITHM: 
!    see formulas in our paper Appendix B
! 
SUBROUTINE CRYSTAL_BENT_ML_CALC_C (oe1,psi_h,i_unit) 

implicit none
type(oeSetup),intent(inout)    ::  oe1
complex(kind=skx),intent(in)   :: psi_h
integer(kind=ski),intent(in)   :: i_unit

real(kind=skr)   :: r_lam0 
real(kind=skr)   :: graze,m1,m2  !,bent_rm,bent_rs
real(kind=skr)   :: gamma_0,gamma_h
real(kind=skr)   :: cry_b,cry_t,cry_a

Real(kind=skr)   :: deltay
Real(kind=skr)   :: polarFactor 

real(kind=skr),dimension(3)   :: vi,vo
real(kind=skr),dimension(6,6) :: s_tensor
Real(kind=skr)      :: VS,A1,A2,A3
Real(kind=skr)      :: csmall,dalpha_over_dt
integer(kind=ski)   :: i_polarFactor

!
! retrieve basic parameters
!
R_LAM0 = TWOPI/oe1%Q_PHOT_SET                    !wavelength in cm

!bragg angle in rad
graze = ASIN(R_LAM0/oe1%crystalData%d_spacing/2.0D0)

m1 = oe1%bent_m1                      ! moment M1
m2 = oe1%bent_m2                      ! moment M2

cry_b = oe1%ass_fac

s_tensor = oe1%crystalElasticity%s

!
! now calculate the dimensionless curvature "c" following formulas of our paper
!
vo = oe1%vout_bragg
vi = oe1%vin_bragg
gamma_0 = vi(3)
gamma_h = vo(3)

! equation B.5 
A1 =-vo(2)*vi(2)*(vo(3)-vi(3))+ &
     vo(2)*vi(3)*(vo(2)-vi(2))+ &
     vo(3)*vi(2)*(vo(2)-vi(2))
A2 = vo(3)*vi(3)*(vo(3)-vi(3))
A3 = vo(3)*vi(3)*(vo(2)-vi(2))

VS = A1*(s_tensor(2,1)*m1+s_tensor(2,2)*m2)+ &
     A2*(s_tensor(3,1)*m1+s_tensor(3,2)*m2)+ &
     A3*(s_tensor(4,1)*m1+s_tensor(4,2)*m2)

! this is d(alphaZ)/dt
dalpha_over_dt = (-2d0*VS/gamma_0)

!
! this is the y-shift between lamellae (different for Bragg and laue)
!
if (oe1%f_refrac.eq.0) then ! bragg
    deltay=2.0D0
else if (oe1%f_refrac.eq.1) then ! laue
    deltay=pi/2.0D0
endif
!
! introduce the direction of the lamella rotation. 
! if concave (R>0) the y of the successive lamellae increases
! if convex (R<0) the y of the successive lamellae reduces
!
deltay = sign(deltay,oe1%bent_rm)

if (i_unit.ge.1) then 
   write(i_unit,*) ' '
   write(i_unit,*) 'A1: ',A1
   write(i_unit,*) 'A2: ',A2
   write(i_unit,*) 'A3: ',A3
   write(i_unit,*) 'VS=-(gamma0/2) d(alpha_zac)/dt: ',VS
   write(i_unit,*) 'd(alpha_zac)/dt: ',dalpha_over_dt
endif

!
! loop over polarization factor
!
do i_polarFactor=0,1
    !polarization factor
    if (i_polarfactor.eq.0) then
      polarFactor = 1.00D0  !s polarized
    else if (i_polarFactor.eq.1) then
      polarFactor = abs(cos(2.00d0*graze))!p polarized
    endif
    
    ! debugging
    if (i_unit.ge.1) then 
      if (i_polarFactor .eq. 0) then
        write(i_unit,*)'  '
        write(i_unit,*)'Sigma-polarization: '
      else 
        write(i_unit,*)'  '
        write(i_unit,*)'Pi-polarization: '
      endif
      write(i_unit,*)'   Polarization factor: ',polarFactor
    endif
    
    ! c following our eq. B.1
    !csmall = dalpha_over_dt*cry_b* & 
    ! R_LAM0*abs(gamma_0)/(polarFactor*cdabs(psi_h)*sqrt(abs(cry_b))*2.0d0*PI)/ &
    !                     (polarFactor*cdabs(psi_h)*sqrt(abs(cry_b)))
    !csmall_b1 = dalpha_over_dt * (cry_b/abs(cry_b)) * abs(gamma_0) * r_lam0/ &
    !                          (2d0*pi*(polarFactor*cdabs(psi_h))**2)
    !
    csmall = dalpha_over_dt*cry_b* & 
                         (0.5d0*oe1%l_ext_s/polarFactor) / &
                         (polarFactor*cdabs(psi_h)*sqrt(abs(cry_b)))
    
    !this is DeltaA
    !TODO:  check sign, removed abs() ?
    cry_a = deltay/abs(csmall)
    
    ! cry_t is the slab thickness, Delta T
    cry_t = cry_a*((oe1%l_ext_s)/polarFactor)
    
    if (i_unit.gt.1) then
        write(i_unit,*) '   Delta(eta): ',deltay
        write(i_unit,*) '   Reduced thickness of a lamella Delta(A): ',cry_a
        write(i_unit,*) '   Thickness of a lamella Delta(T) [cm]: ',cry_t
        write(i_unit,*) '   Number of lamella: ',int(oe1%thickness/cry_t)
        write(i_unit,*) '   c: ',csmall
        write(i_unit,*) '   beta_ml: ',csmall/oe1%l_ext_s
        write(i_unit,*) '    '
        write(i_unit,*) '   D(theta) [eta]: ', csmall*oe1%thickness/oe1%l_ext_s
        write(i_unit,*) '   D(theta) [rad]: ', (csmall*oe1%thickness/oe1%l_ext_s)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)
        !energy resolution: 
        write(i_unit,*) '   DE/E = ', (csmall*oe1%thickness/oe1%l_ext_s)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)/ tan(oe1%graze)
        write(i_unit,*) '   DE [eV] = ', (csmall*oe1%thickness/oe1%l_ext_s)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)/ tan(oe1%graze)*oe1%phot_set
    end if
    !store parameters for sigma-polarization only
    if (i_polarFactor .eq. 0) then 
        oe1%ml_c= csmall
    endif
end do ! i_polarFactor


return
END SUBROUTINE CRYSTAL_BENT_ML_CALC_C



!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!C+++
!C       SUBROUTINE      CRYSTAL_BENT_PP_CALC_BETA
!C
!C       PURPOSE         Computes the beta parameter of the ML model
!C
!C       MODIFIED        Nicolas Perez Bocanegra 2012. 
!C
!C---
!C
!C INPUT PARAMETERS:
!C   OE1:     Type with crystal info
!C   PSI_H:       | Electric susceptibility (proportional to structure factor)
!C
!C OUTPUT PARAMETERS:
!C   sets the parameters:
!C        oe1%bent_G
!C        oe1%bent_beta
!C        oe1%bent_betac
!C
!C ALGORITHM: 
!C   see formulas in our paper Appendix C
!C
SUBROUTINE CRYSTAL_BENT_PP_CALC_BETA ( oe1, PSI_H, I_UNIT) 

implicit none

type(oeSetup),                intent(inout)   ::  oe1
complex(kind=skx),            intent(in)   :: psi_h
integer(kind=ski),            intent(in)   :: i_unit

real(kind=skr),dimension(6,6) :: s_tensor
real(kind=skr),dimension(3)   :: vin1,vout1
real(kind=skr),dimension(3)   :: vtemp,vtemp2
Real(kind=skr)   :: d_spacing,chi,G1,G2,G3
Real(kind=skr)   :: polarFactor,betac, beta 
real(kind=skr)   :: r_lam0,phot,G,m1,m2
real(kind=skr)   :: graze  
real(kind=skr)   :: gamma_0,gamma_h
real(kind=skr)   :: cry_b, pend_len
integer(kind=ski)  :: i_polarFactor

!
! retrieve basic parameters
!
R_LAM0 = TWOPI/oe1%Q_PHOT_SET                    !wavelength in cm
!bragg angle in rad
graze = ASIN(R_LAM0/oe1%crystalData%d_spacing/2.0D0)

m1 = oe1%bent_m1
m2 = oe1%bent_m2

cry_b = oe1%ass_fac

s_tensor = oe1%crystalElasticity%s

!
! Equation for G parameter in Penning-Polder theory and calculation of 
! strain gradient beta
!

!use vin and vout for the exact Bragg directions
vin1 = oe1%vin_bragg
vout1 = oe1%vout_bragg
gamma_0 = vin1(3)
gamma_h = vout1(3)

! Eq. C.6
!G1 = vin1(3)*vout1(3)*bh(3)/sqrt(mod2_bh)
!G2 = vin1(3)*vout1(3)*bh(2)/sqrt(mod2_bh)
!G3 = (vin1(2)*vout1(3)*bh(2) + &
!      vin1(3)*vout1(2)*bh(2) - &
!      vin1(2)*vout1(2)*bh(3) ) /sqrt(mod2_bh)
! Eq. C.6
G1 = vin1(3)*vout1(3)*oe1%bh0(3)
G2 = vin1(3)*vout1(3)*oe1%bh0(2)
G3 = (vin1(2)*vout1(3)*oe1%bh0(2) + &
      vin1(3)*vout1(2)*oe1%bh0(2) - &
      vin1(2)*vout1(2)*oe1%bh0(3) )

! Eq. C.7
G = (1/oe1%CrystalData%d_spacing)*( &
                     G1*( m2*s_tensor(2,3)+m1*s_tensor(1,3) ) +  &
                     G2*( m2*s_tensor(2,4)+m1*s_tensor(1,4) ) +  &
                     G3*( m2*s_tensor(2,2)+m1*s_tensor(1,2) )  )

if (i_unit.ge.1) then 
  write(i_unit,*) '    G1,G2,G3,G: ',G1,G2,G3
  write(i_unit,*) '    G: ',G
endif

!
! this is a loop in polarization factor
!
do i_polarFactor=0,1
    !polarization factor
    if (i_polarFactor.eq.0) then
      polarFactor = 1.00D0  !s polarized
    else if (i_polarFactor.eq.1) then
      polarFactor = abs(cos(2.00d0*graze))!p polarized
    endif

    if (i_unit.ge.1) then 
      if (i_polarFactor .eq. 0) then
        write(i_unit,*)'  '
        write(i_unit,*)'Sigma-polarization: '
      else 
        write(i_unit,*)'  '
        write(i_unit,*)'Pi-polarization: '
      endif
      write(i_unit,*)'   Polarization factor: ',polarFactor
    endif

    !pendellosung length in intensity
    pend_len = r_lam0*abs(gamma_0)/(polarFactor*abs(psi_h))/sqrt(abs(cry_b))
    
    ! this is beta as defined by Xianbo
    beta=G*r_lam0/(polarFactor*abs(psi_h)*sqrt(gamma_0*gamma_h))

    betac = pi/(2.0D0*pend_len) !critical strain gradient
    
    if (i_unit.ge.1) then 
      write(i_unit,*) '    beta [cm^-1]: ',beta
      write(i_unit,*) '    betac [cm^-1]: ',betac
      write(i_unit,*) '    fraction of intensity for new beams: ', &
           exp(-2*PI*abs(betac/beta))
      write(i_unit,*) '    c_pp: ',oe1%l_ext_s*beta
      write(i_unit,*) '    D theta[eta]: ',beta*oe1%thickness
      write(i_unit,*) '    D theta[rad]: ',(beta*oe1%thickness)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)
      write(i_unit,*) '    DE/E = ', (beta*oe1%thickness)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)/ tan(oe1%graze)
      write(i_unit,*) '    DE [eV] = ', (beta*oe1%thickness)*cdabs(psi_h)/sqrt(abs(cry_b))/sin(2.0*oe1%graze)/ tan(oe1%graze) * oe1%PHOT_SET
    endif

    if (i_polarFactor .eq. 0) then 
        oe1%pp_G = G
        !oe1%pp_beta = beta
        !oe1%pp_betac = betac
    endif
    
enddo ! i_polarFactor

RETURN
END SUBROUTINE CRYSTAL_BENT_PP_CALC_BETA

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


!C+++
!C       SUBROUTINE      CRYSTAL_BENT_MOMENTS
!C
!C       PURPOSE         Computes the bending moments from crystal radii
!C
!C       MODIFIED        M Sanchez del Rio
!C
!C---
!C
!C INPUT PARAMETERS:
!C   OE1:     Type with crystal info, in particular:
!C   I_UNIT:  >0 write output in unit=i_unit
!C
!C OUTPUT PARAMETERS:
!C   sets the following parameters
!C       oe1%bent_rm_calc
!C       oe1%bent_rt_calc
!C       oe1%bent_m1
!C       oe1%bent_m2
!C
!C ALGORITHM: 
!C   see formulas in our paper Appendix A
!C
SUBROUTINE CRYSTAL_BENT_MOMENTS ( oe1, I_UNIT) 

implicit none

type(oeSetup),                intent(inout)   ::  oe1
integer(kind=ski),            intent(in)      :: i_unit

real(kind=skr),dimension(6,6) :: s_tensor
real(kind=skr)   :: m1,m2
real(kind=skr)   :: bent_rm,bent_rs
Real(kind=skr)   :: s_along,s_perp

bent_rs=oe1%bent_rs      !sagittal radius in cm
bent_rm=oe1%bent_rm      !meridional/tangential radius in cm

s_tensor = oe1%crystalElasticity%s

!
! If a radius is set to 0 in input, it means that the corresponding moment 
! is zero 
! m1 and m2 are torque over inertia moment
!
if( abs(oe1%bent_rm).le.1d-10 .and. abs(oe1%bent_rs).gt.0) then
  ! SAGITTAL SETUP
  m2=0.0d0
  m1=1.0d0/(oe1%bent_rs*s_tensor(1,1))
  !this is anticlastic radius along MERIDIONAL direction
  bent_rm=1/(m1*s_tensor(2,1)) 
  s_along = s_tensor(1,1)
  s_perp = s_tensor(2,1)
  if (i_unit.ge.1) then 
     write(i_unit,*) 'Sagittal setup: '
     write(i_unit,*) '  Sagittal radius (R1) [cm] : ',oe1%bent_rs
     write(i_unit,*) '  Anticlastic radius (R2) [cm]: ',bent_rm
     write(i_unit,*) '  s(1,1) = s_along: ',s_tensor(1,1)
     write(i_unit,*) '  s(2,1) = s_perp:  ',s_tensor(2,1)
     write(i_unit,*) '  Poison ratio (-s_perp/s_along): ',-s_perp/s_along
     write(i_unit,*) '  M1/I = ',m1
     write(i_unit,*) '  M2/I = ',m2
     write(i_unit,*) 'Using tensor indices (for strain gradient):'
     write(i_unit,*) '  s(2,1) s(3,1) s(4,1): ',s_tensor(2,1),s_tensor(3,1),s_tensor(4,1)
  endif
elseif( abs(oe1%bent_rs).le.1d-10 .and. abs(oe1%bent_rm).gt.0) then
  ! TANGENTIAL SETUP
  m1=0.0d0
  m2=1.0d0/(oe1%bent_rm*s_tensor(2,2))
  !this is anticlastic radius along SAGITTAL direction
  bent_rs=1/(m2*s_tensor(1,2))
  s_along = s_tensor(2,2)
  s_perp = s_tensor(1,2)
  if (i_unit.ge.1) then 
     write(i_unit,*) 'Tangential or meridional setup: '
     write(i_unit,*) '  Tangential or meridional radius (R2) [cm]: ',oe1%bent_rm
     write(i_unit,*) '  Anticlastic radius (R1) [cm]: ',bent_rs
     write(i_unit,*) '  s(2,2) = s_along: ',s_tensor(2,2)
     write(i_unit,*) '  s(1,2) = s_perp:  ',s_tensor(1,2)
     write(i_unit,*) '  Poison ratio (-s_perp/s_along): ',-s_perp/s_along
     write(i_unit,*) '  M1/I = ',m1
     write(i_unit,*) '  M2/I = ',m2
     write(i_unit,*) 'Using tensor indices (for strain gradient):'
     write(i_unit,*) '  s(2,2) s(3,2) s(4,2): ',s_tensor(2,2),s_tensor(3,2),s_tensor(4,2)
  endif
else
  !TANGENTIAL+SAGITTAL SETUP
  !note that bent_rs corresponds to R1 amd bent_rm to R2
  m1=(1/(s_tensor(1,2)*s_tensor(2,1)-s_tensor(1,1)*s_tensor(2,2)))*&
     ((s_tensor(1,2)/oe1%bent_rm)-(s_tensor(2,2)/oe1%bent_rs))
  m2=(1/(s_tensor(1,2)*s_tensor(2,1)-s_tensor(1,1)*s_tensor(2,2)))*&
     ((s_tensor(2,1)/oe1%bent_rs)-(s_tensor(1,1)/oe1%bent_rm))
  !take the average
  s_along = 0.5*(s_tensor(1,1)+s_tensor(2,2))
  s_perp =  0.5*(s_tensor(1,2)+s_tensor(2,1))
  if (i_unit.ge.1) then 
     write(i_unit,*) 'Mixed (Tangential+Sagittal) setup:'
     write(i_unit,*) '  Tangential or meridional radius (R2) [cm]: ',oe1%bent_rm
     write(i_unit,*) '  Sagittal radius (R1) [cm] : ',oe1%bent_rs
     write(i_unit,*) '  s(1,1) s(2,1) s(2,2): ',s_tensor(1,1),s_tensor(2,1),s_tensor(2,2)
     write(i_unit,*) '  s_along = 0.5*(s(1,1)+s(2,2)) = ',s_along
     write(i_unit,*) '  s_perp = 0.5*(s(1,2)+s(2,1)) = ',s_perp
     write(i_unit,*) '  Poison ratio (-s_perp/s_along): ',-s_perp/s_along
     write(i_unit,*) '  M1/I = ',m1
     write(i_unit,*) '  M2/I = ',m2
     write(i_unit,*) 'Using tensor indices (for strain gradient):'
     write(i_unit,*) '  s(2,1) s(3,1) s(4,1): ',s_tensor(2,1),s_tensor(3,1),s_tensor(4,1)
     write(i_unit,*) '  s(2,2) s(3,2) s(4,2): ',s_tensor(2,2),s_tensor(3,2),s_tensor(4,2)
  endif
endif

!-! if (i_unit.ge.1) then 
!-!   write(i_unit,*) '  '
!-!   write(i_unit,*) 'Using tensor indices:'
!-!   write(i_unit,*) '  s(2,1) s(3,1) s(4,1): ',s_tensor(2,1),s_tensor(3,1),s_tensor(4,1)
!-!   write(i_unit,*) '  s(2,2) s(3,2) s(4,2): ',s_tensor(2,2),s_tensor(3,2),s_tensor(4,2)
!-!   write(i_unit,*) 'Rm,Rs [cm]: ',bent_rm,bent_rs
!-! endif

! sets parameters    
    
oe1%bent_rm_calc = bent_rm
oe1%bent_rs_calc = bent_rs
oe1%bent_m1 = m1
oe1%bent_m2 = m2

RETURN
END SUBROUTINE CRYSTAL_BENT_MOMENTS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


END MODULE crystal3
