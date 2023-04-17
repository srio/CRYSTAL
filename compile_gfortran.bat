del *.o *.mod



c:\mingw64\bin\gfortran.exe -static -ffree-line-length-none -O2 -c -cpp -D_COMPILE4WIN shadow_globaldefinitions.f90  
c:\mingw64\bin\gfortran.exe -static -ffree-line-length-none     -c -cpp -D_COMPILE4WIN stringio.f90     
c:\mingw64\bin\gfortran.exe -static -ffree-line-length-none     -c -cpp -D_COMPILE4WIN shadow_math.f90 
c:\mingw64\bin\gfortran.exe -static -ffree-line-length-none     -c -cpp -D_COMPILE4WIN elasticity.f90             
c:\mingw64\bin\gfortran.exe -static -ffree-line-length-none -O2 -c -cpp -D_COMPILE4WIN crystal3.f90
c:\mingw64\bin\gfortran.exe -static -ffree-line-length-none -O2 -c -cpp -D_COMPILE4WIN diff_pat.f90
c:\mingw64\bin\gfortran.exe -static -ffree-line-length-none -O2 -o diff_pat.exe *.o
