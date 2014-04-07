!----
!----
!---- MODULE:  elasticity
!----
!---- crystal elasticity parameters (compliance tensor)
!----
!----
!----

Module elasticity
!---- Use Modules ----!

use shadow_globaldefinitions
use shadow_math, only: dot, cross
use stringio ! in elasticity

!---- Variables ----!
implicit none


!public

!todo move to a common place consistenly with shadow DONE
!real(kind=skr),parameter,public   ::  pi= 3.141592653589793238462643D0
!real(kind=skr),parameter,public   ::  torad= 0.017453292519943295769237D0


type, public, bind(C) :: crystalElasticity
  ! input type: 0:from hkl, 1:from valong,vperp,vnor, 2:from file
  integer(kind=ski)                       :: ielasticity=0  
  ! Poisson ratio (for isotropic crystals)
  real(kind=skr)                          :: poisson=0.0 
  ! tells if elasticity info has been loaded
  !  integer(kind=ski)                       :: loadedflag=0   ! 0:No,1:Yes
  ! crystal 0:Si,1:Ge,2:Diam
  integer(kind=ski)                       :: crystalindex=0 ! 0-2:Si,3:Ge,4:Diam
  ! asymmetry angle in degrees
  real(kind=skr)                          :: alpha=0 
  ! the vectors, 0 means before applying rotation (alpha)
  real(kind=skr),dimension(3)             :: hkl,vnorm,valong,vperp
  real(kind=skr),dimension(3)             :: vnorm0,valong0,vperp0
  ! compliance tensor
  real(kind=skr),dimension(6,6)           :: s 
  ! file name (if ielasticity=3)
  character(len=sklen)                    :: fileElasticity
end type crystalElasticity


!---- Everything is private unless explicitly made public ----!
private 

!---- List of public routines ----!
public :: elasticity_calc_s,elasticity_calc_default
public :: elasticity_prompt,elasticity_calc,elasticity_scan
public :: elasticity_apply_asymmetry,elasticity_report
public :: Rodrigues,get_crystal_s_iso, get_crystal_s

!---- List of private routines ----!

Contains
!
!---- Routines ----!
!
!TODO: move to shadow_math
!C+++++
!C       SUBROUTINE RODRIGUES
!C
!C
!C       PURPOSE Rotates a given vecot by a given angle abot a given axis and
!C               produces a new vector in its place
!C
!C       Note that the angle is in rad. Positive angle is in the "right hand
!C       sense", i.e., in "screw" direction in the direction of axis. 
!C
!C
!C       This is similar to ROTVECTOR in shadow_math.f90, but the rotation
!C       sense os opposite there. 
!C
!C       ALGORTIHM M. Sanchez del Rio, X.Shi, V.Honkimaki & N. Perez Bocanegra 
!C                 2012
!C       http://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
!C
!C
subroutine RODRIGUES (v1,axis,angle,v2)
implicit none
real(kind=skr), dimension(3), intent(in) ::v1,axis
real(kind=skr), intent(in)               :: angle
real(kind=skr), dimension(3), intent(out):: v2
real(kind=skr), dimension(3)             ::vtemp1,axisN
real(kind=skr)                           ::temp2

axisN=axis/sqrt(sum(axis*axis))
call CROSS(axisN,v1,vtemp1)
!print*,'<><>axisN,v1: ',axisN,v1
!print*,'<><>vtemp1: ',vtemp1
call dot(axisN,v1,temp2)
v2=v1*cos(angle)+vtemp1*sin(angle)+axisN*temp2*(1-cos(angle))

return
end subroutine RODRIGUES

!
!---------------------
!

!C+++
!        SUBROUTINE    GET_CRYSTAL_S
!
!        PURPOSE       Obtains crystal elastic compliance tensor from a file
!                      and returns certain crystal parameters required for the
!                      calculation of "c" (ML code) and "G" (PP code)
!
!        CREATED       8 August 2012
!C---

SUBROUTINE GET_CRYSTAL_S(s_tensor)

IMPLICIT none
!type(),intent(in) :: crystal type, Ge, Si
integer                           :: i
real(kind=skr), dimension(6,6), intent(Out) :: s_tensor
s_tensor=0.0
!run code to calculate parameters
!read printed file to transfer parameters into s_tensor
do i=1,2
s_tensor(i,i)=0.592
enddo
s_tensor(3,3)=0.768
do i=4,5
     s_tensor(i,i)=1.26
enddo
s_tensor(6,6)=1.964

s_tensor(1,2)=-0.038
s_tensor(2,1)=-0.038
s_tensor(1,3)=-0.214
s_tensor(3,1)=-0.214
s_tensor(2,3)=-0.214
s_tensor(3,2)=-0.214

RETURN

END SUBROUTINE GET_CRYSTAL_S

!
!---------------------
!
SUBROUTINE GET_CRYSTAL_S_ISO(s_tensor)

IMPLICIT none
!type(),intent(in) :: crystal type, Ge, Si
integer                           :: i
real(kind=skr), dimension(6,6), intent(Out) :: s_tensor
s_tensor=0.0
!run code to calculate parameters
!read printed file to transfer parameters into s_tensor
do i=1,3
s_tensor(i,i)=0.768
enddo

do i=4,6
s_tensor(i,i)=1.26
enddo


s_tensor(1,2)=-0.214
s_tensor(2,1)=-0.214
s_tensor(1,3)=-0.214
s_tensor(3,1)=-0.214
s_tensor(2,3)=-0.214
s_tensor(3,2)=-0.214


RETURN

END SUBROUTINE GET_CRYSTAL_S_ISO

!
!---------------------
!
!
! from a hkl direction calculate default crystal cut:
! vnorm = hkl
! vperp a perpendicular to hkl
! valong = cross product of vnorm,vperp
!
subroutine elasticity_calc_default(hkl,vnorm,valong,vperp)
implicit none
real(kind=skr),dimension(3),intent(in)    :: hkl
real(kind=skr),dimension(3),intent(out)   :: vnorm,valong,vperp

real(kind=skr),dimension(3)    :: vtmp
real(kind=skr)                 :: angle
integer(kind=ski)              :: nzeros,h,k,l


vtmp=0
where(abs(hkl).le.1d-15) 
  vtmp=1
elsewhere
  vtmp=0
endwhere
nzeros = sum(vtmp)

h = hkl(1)
k = hkl(2)
l = hkl(3)

vnorm = hkl

select case (nzeros)
  case(0) 
    !vperp = (/0,-l,k/)
    ! This is Clemens Schulze's choice 
    vperp = (/ -h,-k,(h*h+k*k)/l /)
  case(1)
    if (h.eq.0) then
      vperp = (/1,0,0/)
    else if (k.eq.0) then
      vperp = (/0,1,0/)
    else
      vperp = (/0,0,1/)
    endif
  case(2)
    if (abs(h).gt.1d-6) then
      vperp = (/0,1,0/)
    else if (abs(k).gt.1d-6) then
      vperp = (/0,0,1/)
    else
      vperp = (/1,0,0/)
    endif
  case default
    print*,'elasticity_calc_default: Impossible hkl: ',hkl 
    stop
end select

! rotation of hkl an angle chi to put it along y-axis
! rotation must be counter-screw (minus sign!)
!  memorandum: 
!    alpha = chi + 90
!    alpha = 180 - alphaXOP
!    ->  chi = 90-alphaXOP
! alpha in the code is alphaXOP
!
!angle =  -(90d0)*torad
!call rodrigues(hkl,vperp,angle,valong)
!call cross(vperp,valong,vnorm)
 
call cross(vnorm,vperp,valong)

call check_screw(vperp,valong,vnorm) 

end subroutine elasticity_calc_default

!
!------------------------
!
subroutine elasticity_calc_s(crystalindex,vnorm,valong,vperp,s)


real(kind=skr),dimension(3),intent(in)    :: vnorm,valong,vperp
real(kind=skr),dimension(3)               :: xin,yin,zin
integer(kind=ski),intent(in)              :: crystalindex
real(kind=skr),dimension(6,6),intent(out) :: s
real(kind=skr),dimension(3)               :: x,y,z
real(kind=skr) :: d1,d2
real(kind=skr) :: s11,s12,s44,sc
real(kind=skr) :: s_11,s_12,s_13,s_14,s_15,s_16,s_22,s_23
real(kind=skr) :: s_24,s_25,s_26,s_33,s_34,s_35,s_36,s_44,s_45,s_46,s_55,s_56,s_66


! test orthogonality

call check_screw(vperp,valong,vnorm)
!d1 = dot_product(vnorm,valong)
!d2 = dot_product(valong,vperp)
!IF ((abs(d1).gt.1d-10).or.(abs(d2).gt.1d-10)) then
!  print*,'elasticity_calc_s: Error: inputs vectors are not orthogonal'
!  print*,'   valong: ',valong
!  print*,'   vnorm: ',vnorm
!  print*,'   vperp: ',vperp
!  print*,'   dot_product(vnorm,valong): ',d1
!  print*,'   dot_product(valong,vperp): ',d2
!  stop
!ENDIF

select case (crystalindex) 
  case(0)
     !print*,'elasticity_calc_s: Using silicon'
     s11=0.768            !Si cm**2/10**2 dyn
     s12=-0.214
     s44=1.26
  case(1)
     !print*,'elasticity_calc_s: Using silicon'
     s11=0.768            !Si cm**2/10**2 dyn
     s12=-0.214
     s44=1.26
  case(2)
     !print*,'elasticity_calc_s: Using silicon'
     s11=0.768            !Si cm**2/10**2 dyn
     s12=-0.214
     s44=1.26
  case(3)
     !print*,'elasticity_calc_s: Using germanium'
     s11=0.972
     s12=-0.266
     s44=1.49
  case(4)
     !print*,'elasticity_calc_s: Using diamond'
     s11=0.0949           !C cm**2/10**2 dyn
     s12=0.00978
     s44=0.173
  !case defaut
  !   print*,'elasticity_calc_s: Error: Undefined crystalindex: ',crystalindex
  !   stop
end select

!
! convert from cm^2/(100 dyn) to m^2/N
!
s11 = s11*10
s12 = s12*10
s44 = s44*10

sc=s11-s12-0.5*s44

! Xianbo's reference system
!xin = vnorm
!yin = valong
!zin = vperp

! XOP+SHADOW's reference system
zin = vnorm  
yin = valong
xin = vperp


!normalize
x = xin/sqrt(sum(xin**2))
y = yin/sqrt(sum(yin**2))
z = zin/sqrt(sum(zin**2))
!

s_11=s11+sc*(x(1)**4+x(2)**4+x(3)**4-1)
s_22=s11+sc*(y(1)**4+y(2)**4+y(3)**4-1)
s_33=s11+sc*(z(1)**4+z(2)**4+z(3)**4-1)

s_13=s12+sc*(x(1)**2*z(1)**2 + x(2)**2*z(2)**2 + x(3)**2*z(3)**2)
s_12=s12+sc*(x(1)**2*y(1)**2 + x(2)**2*y(2)**2 + x(3)**2*y(3)**2)
s_23=s12+sc*(y(1)**2*z(1)**2 + y(2)**2*z(2)**2 + y(3)**2*z(3)**2)

s_14=2.*sc*(x(1)**2*y(1)*z(1)+x(2)**2*y(2)*z(2)+x(3)**2*y(3)*z(3))
s_15=2.*sc*(x(1)**2*x(1)*z(1)+x(2)**2*x(2)*z(2)+x(3)**2*x(3)*z(3))
s_16=2.*sc*(x(1)**2*x(1)*y(1)+x(2)**2*x(2)*y(2)+x(3)**2*x(3)*y(3))

s_24=2.*sc*(y(1)**2*y(1)*z(1)+y(2)**2*y(2)*z(2)+y(3)**2*y(3)*z(3))
s_25=2.*sc*(y(1)**2*x(1)*z(1)+y(2)**2*x(2)*z(2)+y(3)**2*x(3)*z(3))
s_26=2.*sc*(y(1)**2*x(1)*y(1)+y(2)**2*x(2)*y(2)+y(3)**2*x(3)*y(3))

s_34=2.*sc*(z(1)**2*y(1)*z(1)+z(2)**2*y(2)*z(2)+z(3)**2*y(3)*z(3))
s_35=2.*sc*(z(1)**2*x(1)*z(1)+z(2)**2*x(2)*z(2)+z(3)**2*x(3)*z(3))
s_36=2.*sc*(z(1)**2*x(1)*y(1)+z(2)**2*x(2)*y(2)+z(3)**2*x(3)*y(3))


s_56=4.*sc*(x(1)**2*y(1)*z(1)+x(2)**2*y(2)*z(2)+x(3)**2*y(3)*z(3))
s_46=4.*sc*(y(1)**2*x(1)*z(1)+y(2)**2*x(2)*z(2)+y(3)**2*x(3)*z(3))
s_45=4.*sc*(z(1)**2*x(1)*y(1)+z(2)**2*x(2)*y(2)+z(3)**2*x(3)*y(3))

s_44=s44+4.*sc*(y(1)**2*z(1)**2+y(2)**2*z(2)**2+y(3)**2*z(3)**2)
s_55=s44+4.*sc*(x(1)**2*z(1)**2+x(2)**2*z(2)**2+x(3)**2*z(3)**2)
s_66=s44+4.*sc*(x(1)**2*y(1)**2+x(2)**2*y(2)**2+x(3)**2*y(3)**2)

s(1, 1)=s_11
s(1, 2)=s_12
s(1, 3)=s_13
s(1, 4)=s_14
s(1, 5)=s_15
s(1, 6)=s_16

s(2, 1)=s_12
s(2, 2)=s_22
s(2, 3)=s_23
s(2, 4)=s_24
s(2, 5)=s_25
s(2, 6)=s_26

s(3, 1)=s_13
s(3, 2)=s_23
s(3, 3)=s_33
s(3, 4)=s_34
s(3, 5)=s_35
s(3, 6)=s_36

s(4, 1)=s_14
s(4, 2)=s_24
s(4, 3)=s_34
s(4, 4)=s_44
s(4, 5)=s_45
s(4, 6)=s_46

s(5, 1)=s_15
s(5, 2)=s_25
s(5, 3)=s_35
s(5, 4)=s_45
s(5, 5)=s_55
s(5, 6)=s_56

s(6, 1)=s_16
s(6, 2)=s_26
s(6, 3)=s_36
s(6, 4)=s_46
s(6, 5)=s_56
s(6, 6)=s_66

end subroutine elasticity_calc_s


!
!---------------------
!

!
! rotate valong and vnorm vectors (atound vperp) to include the 
! affect of the asymmetry angle
!
! alpha is in degrees
!
subroutine elasticity_apply_asymmetry(alpha,vnorm0,valong0,vperp0,&
                                            vnorm ,valong ,vperp)
real(kind=skr),dimension(3), intent(in) ::  valong0,vperp0,vnorm0
real(kind=skr),              intent(in) ::  alpha
real(kind=skr),dimension(3), intent(out)::  valong,vperp,vnorm 
real(kind=skr),dimension(3)             ::  vtmp1
real(kind=skr)                          ::  angle


call check_screw(vperp0,valong0,vnorm0) 

! rotate valong an angle alpha
vperp = vperp0/sqrt(sum(vperp0*vperp0)) !normalized axis (vperp)


! rotation of an angle alpha (as defined in xop) in the 
! screw direction around vperp to put back valong in the 
! y direction
angle = alpha*torad 


valong = valong0/sqrt(sum(valong0*valong0)) !normalized valong
call rodrigues(valong, vperp , angle, vtmp1)
valong = vtmp1/sqrt(sum(vtmp1*vtmp1) )

! vnorm is perpendicular to both vperp and valong
!call cross(vperp,valong,vnorm)
!vnorm = vnorm/sqrt(sum(vnorm*vnorm))

! same for vnorm
vnorm = vnorm0/sqrt(sum(vnorm0*vnorm0) )
call rodrigues(vnorm, vperp , angle, vtmp1)
vnorm = vtmp1/sqrt(sum(vtmp1*vtmp1) )

call check_screw(vperp,valong,vnorm) 
end subroutine elasticity_apply_asymmetry

!
!---------------------
!
subroutine elasticity_report(elas,funit)

implicit none 
type(crystalElasticity),intent(in)     :: elas
integer(kind=ski),intent(in)           :: funit
integer(kind=ski)                      :: i

write(funit,*) ' '
write(funit,'("input type: "I8)')  elas%ielasticity

         select case (elas%ielasticity)
           case(0) 
              write(funit,'("Poisson ratio: "f11.3)')  elas%poisson
           case(1) 
              write(funit,'("crystalindex: "I8)')  elas%crystalindex
              write(funit,'("hkl: "3f8.1)')  elas%hkl
              write(funit,'("alpha [deg]: "f8.1)')  elas%alpha

              write(funit,*) ' '
              write(funit,*) 'Crystallographic directions not affected by asymmetry:'
              write(funit,'("valong0: "3f8.1)')  elas%valong0
              write(funit,'("vnorm0:  "3f8.1)')  elas%vnorm0
              write(funit,'("vperp0:  "3f8.1)')  elas%vperp0
              write(funit,*) ' '
              write(funit,*) 'Crystallographic directions (rotated according asymmetry, normalized):'
              write(funit,'("valong: "3f8.3)')  elas%valong
              write(funit,'("vnorm:  "3f8.3)')  elas%vnorm
              write(funit,'("vperp:  "3f8.3)')  elas%vperp
           case(2) 
              write(funit,'("crystalindex: "I8)')  elas%crystalindex
              write(funit,'("alpha [deg]: "f8.1)')  elas%alpha

              write(funit,*) ' '
              write(funit,*) 'Crystallographic directions not affected by asymmetry:'
              write(funit,'("valong0: "3f8.1)')  elas%valong0
              write(funit,'("vnorm0:  "3f8.1)')  elas%vnorm0
              write(funit,'("vperp0:  "3f8.1)')  elas%vperp0
              write(funit,*) ' '
              write(funit,*) 'Crystallographic directions (rotated according asymmetry, normalized):'
              write(funit,'("valong: "3f8.3)')  elas%valong
              write(funit,'("vnorm:  "3f8.3)')  elas%vnorm
              write(funit,'("vperp:  "3f8.3)')  elas%vperp

           case(3)
              write(funit,*) ' '
              write(funit,'(a)') 'input file with compliance tensor s: '//trim(elas%fileelasticity)
           case default
              print*,'Error: Undefined entry: ',elas%ielasticity
              stop
         end select


write(funit,*) ' '
write(funit,'(a)') 'compliance tensor s: '
do i=1,6
  write(funit,'(6F10.3)') elas%s(i,:)
end do
write(funit,*) ' '
write(funit,'("Poisson ratio meridional -s(1,2)/s(2,2): "f18.5)')  &
    -elas%s(1,2)/elas%s(2,2)
write(funit,'("Poisson ratio sagittal -s(1,2)/s(1,1): "f18.5)')  & 
    -elas%s(1,2)/elas%s(1,1)
write(funit,*) ' '


end subroutine elasticity_report

!
!---------------------
!
!
! fills crystalElasticity inputs from terminal
!
subroutine elasticity_prompt(elas)
        implicit none
        type(crystalElasticity),intent(inout)     :: elas
!
         Write(6,*)'Elasticity info. Obtain compliance tensor (CT) from: '
         Write(6,*)'[0] Poisson ratio (isotropic crystal)'
         Write(6,*)'[1] hkl (for Si/Ge/C)'
         Write(6,*)'[2] crystallographic directions (for Si/Ge/C)'
         Write(6,*)'[3] external file'
         elas%ielasticity =  irint(' <?>')
    
         select case (elas%ielasticity)
           case(0) 
              elas%poisson =  rnumber('Poisson ratio: ')
           case(1) 
              elas%crystalindex =  irint('CrystalIndex: 0,1,2=Si,3=Ge,4=Diamond: ')
              elas%alpha = rnumber('alpha [deg]: ')

              Write(6,*)'Input h k l:'
              Read(*,*) elas%hkl
           case(2) 
              elas%crystalindex =  irint('CrystalIndex: 0,1,2=Si,3=Ge,4=Diamond: ')
              elas%alpha = 0.0
              elas%alpha = rnumber('alpha [deg]: ')
              elas%valong0 = 0.0
              elas%vperp0 = 0.0
              elas%vnorm0 = 0.0
              Write(6,*)'Input valong: '
              Read(*,*) elas%valong0
              Write(6,*)'Input vnorm: '
              Read(*,*) elas%vnorm0
              Write(6,*)'Input vperp: '
              Read(*,*) elas%vperp0
              if ((elas%vperp0(1).eq.0) .and. (elas%vperp0(2).eq.0) &
                                        .and. (elas%vperp0(3).eq.0)) then
                  call cross(elas%valong0,elas%vnorm0,elas%vperp0)
                  print *,'Calculated vperp0: ',elas%vperp0
              endif
              if ((elas%vnorm0(1).eq.0) .and. (elas%vnorm0(2).eq.0) &
                                        .and. (elas%vnorm0(3).eq.0)) then
                  call cross(elas%vperp0,elas%valong0,elas%vnorm0)
                  print *,'Calculated vnorm0: ',elas%vnorm0
              endif
              if ((elas%valong0(1).eq.0) .and. (elas%valong0(2).eq.0) &
                                        .and. (elas%valong0(3).eq.0)) then
                  call cross(elas%vnorm0,elas%vperp0,elas%valong0)
                  print *,'Calculated valong0: ',elas%valong0
              endif
              call check_screw(elas%vperp0,elas%valong0,elas%vnorm0)
           case(3)
              elas%crystalindex=-1
              elas%alpha=0.0
              elas%hkl=0.0
              elas%valong0=0.0
              elas%vnorm0=0.0
              elas%vperp0=0.0
              elas%valong=0.0
              elas%vnorm=0.0
              elas%vperp=0.0
              elas%fileelasticity=rstring('File name containing compliance tensor: ')
           case default
              print*,'Error: Undefined entry: ',elas%ielasticity
              stop
         end select

end subroutine elasticity_prompt

!
!---------------------
!
!
! fills calculated tags in crystalElasticity structure
!
subroutine elasticity_calc(elas)

        type(crystalElasticity),intent(inout)  :: elas
        integer(kind=ski):: i,iscan,alphaN,itmp
        real(kind=skr)   :: alphamin,alphamax,alphastep

         select case (elas%ielasticity)
           case(0) 
              elas%s = 0.0
              elas%s(1,1) = 1.0
              elas%s(2,2) = 1.0
              elas%s(3,3) = 1.0
              elas%s(1,2) = -elas%poisson
              elas%s(2,1) = -elas%poisson
              elas%s(2,3) = -elas%poisson
              elas%s(3,2) = -elas%poisson
           case(1) 
              call elasticity_calc_default(&
                     elas%hkl,&
                     elas%vnorm0, elas%valong0, elas%vperp0)
              call elasticity_apply_asymmetry(elas%alpha, &
                     elas%vnorm0, elas%valong0, elas%vperp0, &
                     elas%vnorm,  elas%valong,  elas%vperp)
              call elasticity_calc_s(& 
                     elas%crystalindex,&
                     elas%vnorm, elas%valong, elas%vperp,&
                     elas%s)
           case(2) 
              call elasticity_apply_asymmetry(elas%alpha, & 
                     elas%vnorm0, elas%valong0, elas%vperp0, &
                     elas%vnorm, elas%valong, elas%vperp )
              call elasticity_calc_s(& 
                     elas%crystalindex,&
                     elas%vnorm, elas%valong, elas%vperp,&
                     elas%s)
           case(3)
              open (20,FILE=elas%fileElasticity,STATUS='OLD')
              read(20,*) elas%s
              close(20)
           case default
              print*,'Error: Undefined entry: ',elas%ielasticity
              stop
         end select

end subroutine elasticity_calc

!
!---------------------
!
!
! makes a scan of the elastic constants as a function of the asymmetry angle
!
subroutine elasticity_scan(elas,alphamin,alphamax,alphan)

        implicit none

        type(crystalElasticity),intent(in)       :: elas
        real(kind=skr),intent(in)                :: alphamin,alphamax
        integer(kind=ski),intent(in)             :: alphaN

        real(kind=skr)   :: alphastep
        integer(kind=ski):: i,iscan,iunit
        type(crystalElasticity)                  :: elas_localcopy


iunit = 6
elas_localcopy = elas
elas_localcopy%alpha  = 0.0

   call elasticity_apply_asymmetry(elas_localcopy%alpha, & 
          elas_localcopy%vnorm0, elas_localcopy%valong0, elas_localcopy%vperp0, &
          elas_localcopy%vnorm, elas_localcopy%valong, elas_localcopy%vperp )
   call elasticity_calc_s(& 
          elas_localcopy%crystalindex,&
          elas_localcopy%vnorm, elas_localcopy%valong, elas_localcopy%vperp,&
          elas_localcopy%s)
   !call elasticity_report(elas_localcopy,iunit)

OPEN   (20,FILE='compliance.spec',STATUS='UNKNOWN')
write(20,'(a)') '#F compliance.spec'
write(20,'(a)') '#S 1 scan alpha '
write(20,'(a)') '#N 24'
write(20,'(a)') '#L alphaX  alpha  chi  s11  s12  s13  s14  s15  s16  s22  s23  s24  s25  s26  s33  s34  s35  s36  s44  s45  s46  s55  s56  s66'
alphastep = (alphamax-alphamin)/(alphan-1)
do i=1,alphan
   elas_localcopy%alpha = alphamin + (i-1)*alphastep
   call elasticity_apply_asymmetry(elas_localcopy%alpha, & 
          elas_localcopy%vnorm0, elas_localcopy%valong0, elas_localcopy%vperp0, &
          elas_localcopy%vnorm, elas_localcopy%valong, elas_localcopy%vperp )
   call elasticity_calc_s(& 
          elas_localcopy%crystalindex,&
          elas_localcopy%vnorm, elas_localcopy%valong, elas_localcopy%vperp,&
          elas_localcopy%s)
   !write(20,'(a)') '#L alpha  s33  s13(23)  s35(34)  s34(35)'
   !write(76,*) elas_localcopy%alpha, elas_localcopy%valong, elas_localcopy%vnorm
   write(20,*) elas_localcopy%alpha, &
               180d0-elas_localcopy%alpha, &
               90d0-elas_localcopy%alpha,&
                   elas_localcopy%s(1,1), &
                   elas_localcopy%s(1,2), &
                   elas_localcopy%s(1,3), &
                   elas_localcopy%s(1,4), &
                   elas_localcopy%s(1,5), &
                   elas_localcopy%s(1,6), &
                   elas_localcopy%s(2,2), &
                   elas_localcopy%s(2,3), &
                   elas_localcopy%s(2,4), &
                   elas_localcopy%s(2,5), &
                   elas_localcopy%s(2,6), &
                   elas_localcopy%s(3,3), &
                   elas_localcopy%s(3,4), &
                   elas_localcopy%s(3,5), &
                   elas_localcopy%s(3,6), &
                   elas_localcopy%s(4,4), &
                   elas_localcopy%s(4,5), &
                   elas_localcopy%s(4,6), &
                   elas_localcopy%s(5,5), &
                   elas_localcopy%s(5,6), &
                   elas_localcopy%s(6,6)
   if (elas_localcopy%alpha .eq. 112) call elasticity_report(elas_localcopy,iunit)
end do
close(20)
write (*,*) 'elasticity_scan: File compliance.spec written to disk.'

end subroutine elasticity_scan
!
!---------------------
!


!
!---------------------
!
!
! from a hkl direction calculate default crystal cut:
! vnorm = hkl
! vperp a perpendicular to hkl
! valong = cross product of vnorm,vperp
!
subroutine check_screw(xx,yy,zz)
implicit none
real(kind=skr),dimension(3),intent(in)    :: xx,yy,zz

real(kind=skr),dimension(3)    :: vtmp1,vtmp2,vtmp3
real(kind=skr)                 :: mtmp1,mtmp2


call cross(xx,yy,vtmp1)
vtmp1=vtmp1/sqrt(sum(vtmp1*vtmp1))
vtmp2=zz/sqrt(sum(zz*zz))
vtmp3 = vtmp1-vtmp2
!print*,'<><><><><><><><><><> DIFF1: ',vtmp3,sum(vtmp3*vtmp3)
if (sum(vtmp3*vtmp3) .gt. 1d-6) stop "CHECK_SCREW: lack of orthogonality..."
 
call cross(yy,zz,vtmp1)
vtmp1=vtmp1/sqrt(sum(vtmp1*vtmp1))
vtmp2=xx/sqrt(sum(xx*xx))
vtmp3 = vtmp1-vtmp2
!print*,'<><><><><><><><><><> DIFF2: ',vtmp1-vtmp2,sum(vtmp1*vtmp1)-sum(vtmp2*vtmp2)
if (sum(vtmp3*vtmp3) .gt. 1d-6) stop "CHECK_SCREW: lack of orthogonality..."
 
call cross(zz,xx,vtmp1)
vtmp1=vtmp1/sqrt(sum(vtmp1*vtmp1))
vtmp2=yy/sqrt(sum(yy*yy))
vtmp3 = vtmp1-vtmp2
!print*,'<><><><><><><><><><> DIFF3: ',vtmp1-vtmp2,sum(vtmp1*vtmp1)-sum(vtmp2*vtmp2)
if (sum(vtmp3*vtmp3) .gt. 1d-6) stop "CHECK_SCREW: lack of orthogonality..."
 


end subroutine check_screw



END MODULE elasticity
