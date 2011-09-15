!PATH = expand_path( '~/IDL:' + !PATH )
define_constants,Rd,Rv,cp,Lv,p0,eps,eps_I,g,pi,Tkel
define_my_constants,text,hours_in_day,minutes_in_hour,secs_in_hour,g_to_kg,hPa_to_Pa


delz = 20.   & ext = 'dz20'
Nz   = 225 

psurf = 97000. ;p 1078
time     = [11.5,15.5,18.,19.0,21.5,24.,26] - 11.5 
SHF      = [-30,90,140,140,100,-10,-10]
LHF      = [5,250,450,500,420,180,0]

surf_temp = 298.
rho_air = psurf/(Rd*surf_temp)
time     = time * secs_in_hour
wt0      = SHF / rho_air / cp
wq0      = LHF / rho_air / Lv
Ntime = N_elements(time)

ug = 10
vg = 0
subs = 0
dqtdx = 0.
dqtdy = 0.

fname = 'ls_flux.inp.ARM_'+ext
openw,1,fname
printf,1,'ARM'
printf,1,'      time      wtsurf      wqsurf    thls      qts       press'
printf,1,'      [s]     [K m/s]    [kg m/s]     [K]      [kg/kg]   [pA]

for t=0,Ntime-1 do begin
    printf,1,time (t),wt0(t),wq0(t),1,1e-10,psurf
endfor


time_ls = [11.5,14.5,17.5,20.5,23.5,26] - 11.5
ddtTemp = [0,0,0,-0.08,-0.16,-0.16] + [-0.125,0,0,0,0,-.1]
ddtqv   = [0.08,0.02,0.04,-.1,-.16,-.3]



time_ls = time_ls * secs_in_hour
Ntime_ls = N_elements(time_ls)
ddtTemp = ddtTemp / secs_in_hour
ddtqv   = ddtqv   / secs_in_hour * g_to_kg

for t=0,Ntime_ls-1 do begin
  printf,1
  printf,1
  printf,1,'large scale forcing terms'
  printf,1,'height     ug    vg      wfls    dqtdx    dqtdy    dqtdt    dthlrad'
  printf,1,'#      ',time_ls(t)

  for k=0,Nz-1 do begin
    zf = (k+0.5) * delz
    if zf le 1000 then begin
       printf,1,format = ('(3f10.4,5e13.4)') ,$
       zf, ug,vg,subs,dqtdx,dqtdy,ddtqv(t),ddtTemp(t)
    endif

    if zf ge 2000 then begin
       printf,1,format = ('(3f10.4,5e13.4)'),$
       zf, ug,vg,subs,dqtdx,dqtdy,0,0
    endif

    if zf gt 1000 and zf lt 2000 then begin
       frac = (2000-zf) / 1000.   ;varies between 0 and 1
       q_tend = ddtqv(t) * frac
       t_tend = ddtTemp(t) * frac
       printf,1,format = ('(3f10.4,5e13.4)'),$
        zf,ug,vg,subs,dqtdx,dqtdy,q_tend,t_tend
    endif

  endfor
endfor


close,1

Ndz40 = 112
text = ''
rec6a=fltarr(6)
rec6b = fltarr(6)
openr,1,'prof.inp.ARM_dz40'
fname='prof.inp.ARM_'+ext
openw,2,fname
fname='lscale.inp.ARM_'+ext
openw,3,fname
h1='input file large scale forcing'
h2='height     ug       vg     wfls    dqtdxls   dqtdyls  dqtdtls        dthlrad'

ls_text = '    10.0000    0.0000    0.0000    0.0000    0.0000    0.0000   0.000' 

printf,3,h1
printf,3,h2
for k=0,1 do begin
 readf,1,text
 printf,2,text
endfor
for k=0,Ndz40 do begin
    readf,1,rec6b
    if k eq 0 then begin
        printf,2,rec6b(0)/2.,rec6b(1:5)  ;10 m 
        rec6a = rec6b
        printf,3,format = ('(f8.1,a)'),rec6b(0)/2.,ls_text
    endif

    if k ge 1 then begin
      printf,2,rec6a + (rec6b-rec6a) * 0.25
      printf,2,rec6a + (rec6b-rec6a) * 0.75

      printf,3,format = ('(f8.1,a)'),rec6a(0) + (rec6b(0)-rec6a(0)) * 0.25,ls_text
      printf,3,format = ('(f8.1,a)'),rec6a(0) + (rec6b(0)-rec6a(0)) * 0.75,ls_text
      rec6a = rec6b
    endif

    ;if k eq Ndz40 then begin
    ;  printf,2,rec6b
    ;  printf,3,format = ('(f8.1,a)'),rec6b(0),ls_text
    ;endif
    
endfor
close,1
close,2
close,3



END

