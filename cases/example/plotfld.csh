#! /bin/csh

ncl plotfld.ncl 'fname=(/"tmser.001.nc"/)' \
'ffld=(/"cfrac","zi","H","LE"/)' 'foname="time"' ncols=1 npages=2 tbeg=0. tend=3600.

ncl plotfld.ncl 'fname=(/"profiles.001.nc"/)' \
'ffld=(/"u","v","thl","qt","ql","thv","wtvt"/)' 'foname="profiles"' ncols=2 npages=2 tbeg=0. tend=3600. ymax=1200. ymin=00.

