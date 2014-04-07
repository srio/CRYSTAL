PROGRAM compliance2

use shadow_globaldefinitions
use stringio
!use shadow_math, only: dot, rotvector, norm , cross
use shadow_math, only: cross
use elasticity

        implicit none

        integer(kind=ski):: i,iscan,alphaN,itmp
        real(kind=skr)   :: alphamin,alphamax
        type(crystalElasticity)       :: elas

!        real(kind=skr)   :: angle
!        real(kind=skr),dimension(3)   :: a,v,v2
!
!100 continue
!Write(6,*)'Input v: '
!Read(*,*) v
!Write(6,*)'Input axis: '
!Read(*,*) a
!Write(6,*)'angle: '
!Read(*,*) angle
!angle=angle*pi/180
!
!call rodrigues(v,a,angle,v2)
!Write(6,*)'result: ',v2
!goto 100

call elasticity_prompt(elas)
call elasticity_calc(elas)

itmp = 6
call elasticity_report(elas,6)

! 
! !C
OPEN   (20,FILE='compliance.dat',STATUS='UNKNOWN')
do i=1,6
  write(20,'(6F10.3)') elas%s(i,:)
end do
close(20)
write (*,*) '>> File compliance.dat written to disk.'



iscan = 0
if ((elas%ielasticity .eq. 1) .or. (elas%ielasticity .eq. 2))  iscan = irint('Perform alpha scan? ')
if (iscan.eq.1) then
   alphamin = rnumber('  From (alphamin, deg): ')
   alphamax = rnumber('    To (alphamax, deg): ')
   alphaN = irint('    number of points: ')
   call elasticity_scan(elas,alphamin,alphamax,alphaN)
endif


! 
END PROGRAM compliance2
