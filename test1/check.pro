dir = './data'
nnnn=64
num =400




x = findgen(nnnn)/(nnnn)
y = findgen(nnnn)/(nnnn)
z = findgen(nnnn)/(nnnn)

bx = fltarr(nnnn,nnnn,nnnn)
by = fltarr(nnnn,nnnn,nnnn)
bz = fltarr(nnnn,nnnn,nnnn)
kk = 2*!pi
lambda = 1.0
ll = sqrt(kk^2 - lambda^2)

for i=0,nnnn-1 do begin
for j=0,nnnn-1 do begin
for k=0,nnnn-1 do begin
    bx[i,j,k] = ll*sin(kk*x[i])*cosh(ll*(1.-z[k]))
    by[i,j,k] = lambda*sin(kk*x[i])*sinh(ll*(1.-z[k]))
    bz[i,j,k] = kk*cos(kk*x[i])*sinh(ll*(1.-z[k]))
endfor
endfor
endfor

bx = bx/(max(abs(bz[*,*,0])))
by = by/(max(abs(bz[*,*,0])))
bz = bz/(max(abs(bz[*,*,0])))

b30=dblarr(nnnn,nnnn,nnnn,3)


openr,lun,dir+'/B_0002_0032.dat',/get_lun
readu,lun,b30
free_lun,lun


posi1=[0.03,0.05,0.33,0.95]
posi2=[0.35,0.05,0.65,0.95]
posi3=[0.67,0.05,0.97,0.95]

bx30 = reform(b30[*,*,*,0])
by30 = reform(b30[*,*,*,1])
bz30 = reform(b30[*,*,*,2])

wdef,0,1500,500
cgplot,bx[*,*,*],bx30[*,*,*],psym=1,xr=[-1,1],yr=[-1,1],position=posi1
;cgoplot,[-10,10],[-10,10],linestyle=1,color='red'
cgplot,by[*,*,*],by30[*,*,*],psym=1,xr=[-1,1],yr=[-1,1],position=posi2,/noer
;cgoplot,[-10,10],[-10,10],linestyle=1,color='red'
cgplot,bz[*,*,*],bz30[*,*,*],psym=1,xr=[-1,1],yr=[-1,1],position=posi3,/noer
;cgoplot,[-10,10],[-10,10],linestyle=1,color='red'


;cgplot,bx,transpose(by100,[1,0,2]),psym=1,xr=[-1,1],yr=[-1,1],position=posi1
;cgplot,by,transpose(-1*bx100,[1,0,2]),psym=1,xr=[-1,1],yr=[-1,1],position=posi2,/noer
;cgplot,bz,transpose(bz100,[1,0,2]),psym=1,xr=[-1,1],yr=[-1,1],position=posi3,/noer


;jx = finite_diff(Bz100,2) - finite_diff(By100,3)
;jy = finite_diff(Bx100,3) - finite_diff(Bz100,1)
;jz = current_jz(Bx100,By100,Bz100)

;vec_2vtk,bx,by,bz,'bxyz_real.vtk'
;vec_2vtk,bx100,by100,bz100,'bxyz100.vtk'
;vec_2vtk,jx*100,jy*100,jz*100,'jxyz.vtk'
end

