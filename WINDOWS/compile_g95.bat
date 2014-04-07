del *.o *.mod
set PATH = C:\g95\bin          
set G95_LIBRARY_PATH=C:\g95\lib

g95 -static -ffree-line-length-huge -O2 -fomit-frame-pointer -c shadow_globaldefinitions.F90
g95 -static -ffree-line-length-huge -O2 -fomit-frame-pointer -c stringio.F90                  
g95 -static -ffree-line-length-huge -O2 -fomit-frame-pointer -c shadow_math.F90 
g95 -static -ffree-line-length-huge -O2 -fomit-frame-pointer -c elasticity.F90     
g95 -static -ffree-line-length-huge -O2 -fomit-frame-pointer -c crystal3.F90         
g95 -static -ffree-line-length-huge -O2 -fomit-frame-pointer -c diff_pat.F90
g95 -static -ffree-line-length-huge -O2 -fomit-frame-pointer -o diff_pat.exe *.o