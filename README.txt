
New version of xcrystal and its binary module diff_pat that
integrates now bent crystal models (ML multilamellar and PP Penning-Polder). 

Now xcrystal_bent becomes obsolete (waiting for the TT module...)

In addition, a new main proglam "compliance" is provided to calculate
compliance tensors. It is not interfaced to XOP right now.

srio@esrf.eu    2014-03-24



In Mac: 

make -f makefile.mac clean all

In Linux: 

make -f makefile.unix clean all
