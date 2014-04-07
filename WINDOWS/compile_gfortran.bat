
set PATH=c:/MinGW;%PATH%

del *.o *.mod
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -c shadow_globaldefinitions.F90
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -c stringio.F90   
# -O2 give problems here, why?               
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none     -c shadow_math.F90 
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -c elasticity.F90     
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -c crystal3.F90         
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -c diff_pat.F90
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -o diff_pat.exe *.o

del diff_pat.o
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -c compliance.F90
c:\mingw\bin\gfortran.exe -static -ffree-line-length-none -O2 -o compliance.exe *.o

#
#install
# 

copy diff_pat.exe ..\..\..\bin.x86\diff_pat.exe
copy compliance.exe ..\..\..\bin.x86\compliance.exe

#
#clean 
#

del *.o *.mod