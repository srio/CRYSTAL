!C
!C This program calculates the diffraction pattern of a crystal.
!C

!C
!C MSR/94/03/25 first version hacked from Shadow code
!C MSR/97/10/07 Adapts for xop 1.9 (embedded version)
!C MSR/00/06/13 Adapts for xop 1.0
!C        renamed from diff_pattern to diff_pat
!C MSR/08/09/05 Updated TOCM,TOANGS
!C Changed by Eric.Lebigot@normalesup.org
!C MSR/09/04/14 Updated TOCM,TOANGS everywhere. 
!C srio@esrf.eu 2011-08-18 Ported to f90. Lots of modifications. 
!C
!C 2012-12-07 srio@esrf.eu
!C Added Bent crystal models. Work done with Nicolas Perez-Bocanegra
!C Xianbo Shi and Veijo Honkimaki
!C 2014-01-22 srio@esrf.eu New cleaning.
!C

PROGRAM diff_pat

use shadow_globaldefinitions
use stringio
use shadow_math, only: dot, norm , cross !, rotvector
use elasticity
use crystal3

implicit none

integer(kind=ski):: scan_mode,scanunit,scanpoints,i
real(kind=skr)   :: phot,q_phot
real(kind=skr)   :: r_s, r_p,phase_s, phase_p, depht_mfp_s
real(kind=skr)   :: theta_b, depht_mfp_p
real(kind=skr)   :: scanmin,scanmax,scanstep,theta_energy,graze,graze0
real(kind=skr)   :: scanval,scanval2
real(kind=skr)   :: sin_ref,theta_b_ref,scanval_ref,graze_ref
real(kind=skr)   :: tmp
real(kind=skr)   :: a_braggDeg
Real(kind=skr)   :: inp_flag,bent_phot,bent_geo,bent_rm,bent_rs, bent_method
Real(kind=skr)   :: bent_alph, bent_thick, bent_ymax, bent_ymin, bent_dy0 
Integer          :: ierror,itmp
character(len=sklen) ::ans
!vectors
real(kind=skr),dimension(3)   ::  vvin, vvout, bh, bh0, vnor, axis
real(kind=skr),dimension(3)   ::  vin_bragg,vin_bragg_uncorr
real(kind=skr),dimension(3)   ::  vpar, vin_bragg_energy, vtmp

type(oeSetup)    :: oe1

integer(kind=ski):: iplotangles=0 !plot angles with gle
!debug variables
real(kind=skr)   :: gamma_0,gamma_h
integer(kind=ski):: i_debug= 0 !11 !debug: set the flag to the fort unit

!C
!C Version number also in another place!!!!!!
!C
write(6,*) ' '
write(6,*) ' '
write(6,*) ' '
write(6,*) '*****************************************************'
write(6,*) '*           DIFF_PAT    v1.8  (24 Mar 2014)         *'
write(6,*) '*Calculation of a single crystal diffraction profile*'
write(6,*) '*                                                   *'
write(6,*) '*****************************************************'
write(6,*) ' '
oe1%file_refl = rstring(' Name of file with crystal data (from BRAGG): ')

!C
!C Inquires about perfect or mosaic crystal
!C
write (6,*) 'What kind of crystal you want to use ?: '
write (6,*) '        [0] Perfect crystal '
write (6,*) '        [1] Mosaic crystal '
Write (6,*) '        [2] Bent Multilamellar (ML)'
Write (6,*) '        [3] Bent Penning-Polder (PP)'
!        Write (6,*) '        [4] Bent Takagi-Taupin'
oe1%f_mosaic = irint(' <?> ')

!C
!C Inquires about Reflection (Bragg) or Transmission (Laue) geometry.
!C If the crystal is perfect we also allow to calculate the transmitted
!C (not diffracted beams.
!C
write (6,*) 'What do you want to calculate ?: '
if (oe1%f_mosaic.eq.0) then
    write (6,*) '[0] Diffracted beam in Reflection (Bragg) geometry'
    write (6,*) '[1] Diffracted beam in Transmission (Laue) geometry'
    write (6,*) '[2] Transmitted beam in Bragg case'
    write (6,*) '[3] Transmitted beam in Laue case'
else if (oe1%f_mosaic.gt.0) then
    write (6,*) '[0] Diffracted beam in Reflection (Bragg) geometry'
    write (6,*) '[1] Diffracted beam in Transmission (Laue) geometry'
endif
oe1%f_refrac = irint(' <?> ')

!C
!C Inquires about mosaic crystal values.
!C
if (oe1%f_mosaic.eq.1) then
    oe1%spread_mos  = RNUMBER('mosaic angle spread (FWHM) [deg] ? ')
    oe1%spread_mos = torad*oe1%spread_mos/2.35
else
    oe1%spread_mos   = 0.0d0
end if
!C
!C Inquires crystal Thickness. In the Perfect crystal case, we allow
!C to use the thick crystal approximation. In such case, the user must
!C input any negative value.
!C
WRITE(6,*) 'Input the thickness of the crystal [cm] '
if (oe1%f_mosaic.eq.0) then
    write(6,*) ' [Type a negative value for using the thick crystal approximation] '
endif
oe1%thickness = RNUMBER(' <?> ')

!C
!C Inquires about asymmetrical diffraction
!C
oe1%f_bragg_a = 0
IF (oe1%F_MOSAIC.NE.1) THEN
    WRITE(6,*) 'Asymmetric cut angle (deg) between face and bragg planes (CW)= '
    READ(*,*) A_BRAGGDEG
    if (a_braggDeg.eq.0.0) then
        oe1%f_bragg_a = 0
    else
        oe1%f_bragg_a = 1
    endif
ELSE  ! mosaic case
    if (oe1%f_refrac.EQ.1) a_braggDeg = 90.0
    if (oe1%f_refrac.EQ.0) a_braggDeg = 0.0
ENDIF
oe1%a_bragg = a_braggDeg*torad


if ((oe1%f_mosaic.eq.3).and.(oe1%f_refrac.eq.0)  ) then
    print*,'DIFF_PATT: Error in Bragg Geometry: Penning Polder method only for Laue'
    print*,'           ** aborting run **'
    stop
endif

if ((oe1%f_mosaic.eq.3).and.(oe1%f_refrac.eq.2)  ) then
    print*,'DIFF_PATT: Error in Bragg Geometry: Penning Polder method only for Laue'
    print*,'           ** aborting run **'
    stop
endif


!C
!C Inquires about the scanning variable
!C
write (6,*) 'Please select scanning variable: '
write (6,*) ' [1] Incident/Reflected angle [absolute] '
write (6,*) ' [2] Incident/Reflected angle minus theta Bragg corrected'
write (6,*) ' [3] Incident/Reflected angle minus theta Bragg'
write (6,*) ' [4] Photon energy '
write (6,*) ' [5] y variable [Zachariasen] '

!C
!C Call CRYSTAL to read the file with the crystal data (kwhat<0)
!C
call crystal(-1, oe1, &   ! we only need thes inputs, the 
                            ! following ones are not used here.
            q_phot, vvin, vvout, bh, vnor, &
            r_s, r_p,phase_s, phase_p, depht_mfp_s, depht_mfp_p) 
            ! , theta_b, ssr, spr)

! alternatively we can use:
!Call crystal_loadCrystalData (oe1%FILE_REFL,oe1%crystalData)
!C
!C Define scan mode (Theta or Energy)
!C
scan_mode = irint('<?> ')
if (scan_mode.lt.4) then
    PHOT    = RNUMBER ('... at what energy (eV) ? ')
    Q_PHOT  = PHOT*TWOPI/TOCM       ! 2 pi / lambda
else if (scan_mode.eq.4) then
    theta_energy = RNUMBER ('... at what grazing angle [deg]: ')
    IF (theta_energy.LE.0) THEN
        PHOT = -1.0D0*theta_energy
        Q_PHOT  = PHOT*TWOPI/TOCM       ! 2 pi / lambda
        theta_energy = 0
    ELSE
        theta_energy = theta_energy*torad
        PHOT = 2.0D0*oe1%crystalData%D_SPACING*SIN(theta_energy) ! lambda (cm)
        PHOT = TOCM/PHOT                  ! energy [eV]
        Q_PHOT  = PHOT*TWOPI/TOCM       ! 2 pi / lambda
    ENDIF
else if (scan_mode.eq.5) then
    PHOT    = RNUMBER ('... at what energy (eV) ? ')
    Q_PHOT  = PHOT*TWOPI/TOCM       ! 2 pi / lambda
endif

!C
!C Call CRYSTAL to calculate useful parameters
!C

OPEN    (23,FILE='diff_pat.par',STATUS='UNKNOWN')
REWIND (23)

! call crystal in set mode (kwhat=0)
call crystal(0, oe1, q_phot, vvin, vvout, bh, vnor, & 
             r_s, r_p,phase_s, phase_p, depht_mfp_s, depht_mfp_p) 
             !, theta_b, ssr, spr)
! moved ahead to include elasticity 
!! call crystal in info mode (kwhat=1)
!   call crystal  (1, oe1, q_phot, vvin, vvout, bh, vnor, & 
!r_s, r_p,phase_s, phase_p, depht_mfp_s, depht_mfp_p, theta_b, &
!ssr, spr)

!
! check for impossible geometry
!
if ((oe1%f_refrac.eq.0).or.(oe1%f_refrac.eq.2)  ) then
    if (abs(oe1%a_bragg).gt.abs(oe1%graze)) then
        print*,'DIFF_PATT: WARNING in Bragg Geometry: Bragg angle is usually larger than asymmetry angle'
        print*,'           Bragg angle [deg]: ',oe1%graze*todeg
        print*,'           Asymmetry angle [deg]: ',oe1%a_bragg*todeg
        !print*,'           ** aborting run **'
        !stop
    endif
endif
if ((oe1%vin_bragg_uncorr(3).gt.0)) then
    print*,'DIFF_PATT: WARNING in Geometry: incident beam is pointing outwards the cystal surface (0,0,1)'
endif
!!call plotgle(oe1%vin_bragg_uncorr,oe1%vout_bragg_uncorr,oe1%bh0,oe1,iplotangles)


graze0 = asin(pi/q_phot/oe1%crystalData%d_spacing) 
graze     = graze0 + oe1%a_bragg
graze_ref = graze0 - oe1%a_bragg

!C
!C Units for the scanning variable
!C
if (scan_mode.lt.4) then 
    write(6,*) 'Select units for the scanning variable (angle):'
    write(6,*) '[0] radians'
    write(6,*) '[1] microradians'
    write(6,*) '[2] degrees'
    write(6,*) '[3] arcsec'
    scanunit = irint(' <?>')
endif

write(6,*) 'Input scanning variable limits in the chosen units.'
scanmin = rnumber('Minimum ? ')
scanmax = rnumber('Maximum ? ')
scanpoints = irint('Number of points ? ')
scanstep = (scanmax - scanmin) / (scanpoints - 1)

If(oe1%f_mosaic .gt. 1)then
    Write(6,*)'Input sagittal bending radius [cm].'
    READ(*,*) oe1%bent_rs
    Write(6,*)'Input meridional radius [cm].'
    READ(*,*) oe1%bent_rm
    !Write(6,*)'Input Poisson Ratio.'
    !READ(*,*) oe1%kstar
    !todo: remove, bent_pol is not used.
    !Write(6,*)'Polarisation factor:'
    !Write(6,*)'[0] s'
    !Write(6,*)'[1] p'
    !Write(6,*)'[2] unpolarised'
    tmp = irint(' <?>')
    !
    ! elasticity inputs 
    !
    call elasticity_prompt(oe1%crystalElasticity)
ENdIF

!
! elasticity calculations
!
If(oe1%f_mosaic .gt. 1)then
    call elasticity_calc(oe1%crystalElasticity) 
    itmp = 6
    call elasticity_report(oe1%crystalElasticity,itmp)
endif

! call crystal in info mode (kwhat=1)
WRITE(*,*) ' '
WRITE(*,*) 'So far, we are working with:'
call crystal(1, oe1, q_phot, vvin, vvout, bh, vnor, & 
             r_s, r_p,phase_s, phase_p, depht_mfp_s, depht_mfp_p) 
             !, theta_b, ssr, spr)

! make graphical sketches
call plotgle(oe1%vin_bragg_uncorr,oe1%vout_bragg_uncorr,oe1%bh0,oe1,iplotangles)
call plotidl(oe1%vin_bragg_uncorr,oe1%vout_bragg_uncorr,oe1%bh0)

!

!todo: rm these variables, make direct call
vnor = oe1%vnor
vpar = oe1%vpar
axis = oe1%axis
bh = oe1%bh
bh0 = oe1%bh0
vin_bragg = oe1%vin_bragg
vin_bragg_uncorr = oe1%vin_bragg_uncorr
vin_bragg_energy = oe1%vin_bragg_energy
vvout = oe1%vout_bragg

!------------------------------------------------------------------------------

!write(99,*) 'calling crystal with exact Bragg angle: '
!write(99,*) 'vin_bragg: ',vin_bragg
!write(99,*) 'vvout: ',vvout
call crystal  (2, oe1, q_phot, vin_bragg, vvout, bh, vnor, &
     r_s, r_p,phase_s, phase_p, depht_mfp_s, depht_mfp_p) 
     !, theta_b, ssr, spr)
!Write(99,*) 'reflectivities at y=0: ', r_s,r_p

!C 
!C write output file
!C
OPEN   (20,FILE='diff_pat.dat',STATUS='UNKNOWN')
REWIND (20)
write(20,'(a)') '#F diff_pat.dat'
write(20,'(a)') '#S 1 diff_pat run'
write(20,'(a)') '#C results of diff_pat run'
write(20,'(a)') '#N 7'

if (scan_mode.EQ.4) then 
          write(20,'(a)') '#L  E[eV]  Lambda[A]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized-profile  s-polarized-profile'
else
    if (scan_mode.eq.1) then
        if (scanunit.eq.1) then
            write(20,'(a)') '#L  Theta{in} [microrad]  Theta{out} [microrad]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
        else if (scanunit.eq.2) then
            write(20,'(a)') '#L  Theta{in} [deg]  Theta{out} [deg]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
        else if (scanunit.eq.3) then
            write(20,'(a)') '#L  Theta{in} [arcsec]  Theta{out} [arcsec]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
        endif
    else if (scan_mode.eq.2) then
        if (scanunit.eq.1) then
            write(20,'(a)') '#L  Th-ThBc{in} [microrad]  Th-ThBc{out} [microrad]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
        else if (scanunit.eq.2) then
            write(20,'(a)') '#L  Th-ThBc{in} [deg]  Th-ThBc{out} [deg]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
        else if (scanunit.eq.3) then
            write(20,'(a)') '#L  Th-ThBc{in} [arcsec]  Th-ThBc{out} [arcsec]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
         endif
    else if (scan_mode.eq.3) then
         if (scanunit.eq.1) then
             write(20,'(a)') '#L  Th-ThB{in} [microrad]  Th-ThB{out} [microrad]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
         else if (scanunit.eq.2) then
             write(20,'(a)') '#L  Th-ThB{in} [deg]  Th-ThB{out} [deg]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
         else if (scanunit.eq.3) then
             write(20,'(a)') '#L  Th-ThB{in} [arcsec]  Th-ThB{out} [arcsec]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
         endif
    else if (scan_mode.eq.5) then
             write(20,'(a)') '#L  y  y  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized'
    endif
endif

!C
!C
!C
do i=1,scanpoints !C main loop over the scanning points ===== starts ====
    scanval = scanmin + scanstep*(i-1)
    scanval2 = scanval

    !C
    !C convert user unit to shadow units (rads, eV)
    !C 
    If(scan_mode.lt.4) then
        if (scanunit.eq.1) then
          scanval2 = scanval*1D-6
        else if (scanunit.eq.2) then
          scanval2 = scanval*torad
        else if (scanunit.eq.3) then
          scanval2 = scanval*torad/3600D0
        else
          scanval2 = scanval
        end if
    end IF

    select case (scan_mode)    !scanmode is case type
        ! minus sign in angle is to perform cw rotation when 
        ! scanval increses
        case (1)  !  angle scan absolute
            !call rotvector(vnor, axis, scanval2,vvin)
            call rodrigues(vnor, axis, -scanval2,vvin)
        Case(2) ! angle scan (theta-theta_B_corr)
            !call rotvector(vin_bragg, axis, scanval2, vvin)
            call rodrigues(vin_bragg, axis, -scanval2, vvin)
        case(3) ! angle scan (theta-thetaBragg)
            !call rotvector(vin_bragg_uncorr, axis, scanval2, vvin)
            call rodrigues(vin_bragg_uncorr, axis, -scanval2, vvin)
        case(4) ! energy scan
            vvin = vin_bragg_energy
            q_phot  = scanval2*twopi/tocm
        case(5) ! y scan
            if ((oe1%f_refrac.eq.0).or.(oe1%f_refrac.eq.2)) then ! bragg
                !TODO: check this minus sign!!!! DONE
                call rodrigues(vin_bragg,axis,-scanval2*abs(oe1%ssr),vvin)
            else ! laue
                call rodrigues(vin_bragg,axis,scanval2*abs(oe1%ssr),vvin)
            endif
    end select
    !C
    !C calculate output direction 
    !C
    call scat_direction(vvin,bh,vnor,q_phot,vvout)
    call dot (vvout,vnor,sin_ref)
    call dot (vvin,vnor,tmp)

    !C
    !C calculate reflectivity
    !C
    
    !C
    !srio debugging
    !        write(99,*) ' ' 
    !        write(99,*) '---- ' 
    !        write(99,*) 'ssr: ',ssr
    if (i_debug .ge. 1) then
        Call dot (vnor,vvin,gamma_0)
        Call dot (vnor,vvout,gamma_h)
        write(11,*) ' '
        write(11,*) '-----'
        !note the change of sign because Zachariasen normal is pointing
        !into the crystal, and we point outside.
        write(i_debug,*) 'Zac gamma_0,gamma_h (diff_pat): ',-gamma_0,-gamma_h
        write(11,*) 'scanval: ',scanval
    endif
    call crystal(2, oe1, q_phot, vvin, vvout, bh, vnor, &
                 r_s, r_p,phase_s, phase_p, depht_mfp_s, depht_mfp_p) 
                 !, theta_b, ssr, spr)

    !C
    !C if asymmetric is selected, calculate and write down
    !C the output angle
    !C
    if (oe1%f_bragg_a.EQ.1) then
        call dot (vvout,vnor,sin_ref)
  
        if (scan_mode.EQ.1) then
            scanval_ref = asin(sin_ref)
        else if (scan_mode.EQ.2) then
            scanval_ref =  - theta_b_ref + asin(sin_ref)
        else if (scan_mode.EQ.3) then
            scanval_ref =  - graze_ref + asin(sin_ref)
        else if (scan_mode.EQ.4) then
        !C scanval_ref = scanval
        else if (scan_mode.EQ.5) then
           scanval_ref = scanval
        endif
  
        if (scanunit.eq.1) then
            scanval_ref = scanval_ref/1D-6
        else if (scanunit.eq.2) then
            scanval_ref = scanval_ref/torad
        else if (scanunit.eq.3) then
            scanval_ref = scanval_ref/torad*3600D0
        else
            scanval_ref = scanval_ref
        end if
    else
        scanval_ref = scanval
    endif  
    if (scan_mode.EQ.4) scanval_ref = TOANGS/scanval
    WRITE (20,'(7(g15.8,1x))') scanval,scanval_ref, &
           phase_p,phase_s, r_p*r_s*sin(phase_p-phase_s), &
           r_p**2,r_s**2

enddo !C main loop over the scanning points ===== ends ====

!C
!C close files and exit
!C
close(23)
close(20)
write (*,*) '>> '
write (*,*) '>> Files diff_pat.par (parameters) and '
if (scan_mode.EQ.4) then 
    write (*,*) '>> diff_pat.dat (E [eV],Lambda [A],...,p-pol,s-pol)'
else
    write (*,*) '>> diff_pat.dat (Scan in,Scan out,...,p-pol,s-pol)'
endif
write (*,*) '>> written  to disk.'
write (*,*) '>> '


END PROGRAM diff_pat
