
function xderiv2,arr,DX=dx
  ss=size(arr)
  NX=ss[1] & NY=ss[2]
  ; Assumption of a unit size grid in the x-direction, and dx=dy
  if (arg_present(dx) eq 0) then begin
    dx=1.0/(NX-1.0)
  endif  
  hidx=0.5/dx 
  ; Points in the volume
  darr=hidx*(shift(arr,-1,0)-shift(arr,+1,0))
  ; Points at left hand edge
  darr[0,*]=hidx*(-3.*arr[0,*]+4.*arr[1,*]-arr[2,*])
  ; Points at right hand edge
  darr[NX-1,*]=hidx*(3.*arr[NX-1,*]-4.*arr[NX-2,*]+arr[NX-3,*])
  return,darr
end

; yderiv2.pro
;
; Calculate the derivative in the y-direction of an NX x NY scalar field. 
; Three point differencing (error is step size squared)
function yderiv2,arr,DY=dy
  ss=size(arr)
  NX=ss[1] & NY=ss[2] 
  ; Assumption of a unit size grid in the x-direction, and dx=dy
  if (arg_present(dy) eq 0) then begin
    dx=1.0/(NX-1.0)
    dy=dx
  endif  
  hidy=0.5/dy 
  ; Points in the volume
  darr=hidy*(shift(arr,0,-1)-shift(arr,0,+1))
  ; Points at left hand edge
  darr[*,0]=hidy*(-3.*arr[*,0]+4.*arr[*,1]-arr[*,2])
  ; Points at right hand edge
  darr[*,NY-1]=hidy*(3.*arr[*,NY-1]-4.*arr[*,NY-2]+arr[*,NY-3])
  return,darr
end


; Storruck linear force-free field
dim = 64

x = ( (dblarr(dim)*0+1)##findgen(dim) )/(dim)
y = ( findgen(dim)##(dblarr(dim)*0+1) )/(dim)

kk = 2*!pi
lambda = 1.0
ll = sqrt(kk^2 - lambda^2)

bx = ll*sin(kk*x)*cosh(ll*(1.))
by = lambda*sin(kk*x)*sinh(ll*(1.))
bz = kk*cos(kk*x)*sinh(ll*(1.))

dxBy = xderiv2(by)/dim
dyBx = yderiv2(bx)/dim

jz = dxBy - dyBx
alpha = (dim-1)*jz/bz
alpha = alpha*0+median(alpha)

bx = bx/(max(abs(bz)))
by = by/(max(abs(bz)))
bz = bz/(max(abs(bz)))

OPENW, 1, './data/bz0.dat'
WRITEU, 1, double(bz)
close,1

OPENW, 1, './data/alpha0.dat'
WRITEU, 1,double(alpha)
close,1

OPENW, 1, './data/sig0.dat'
WRITEU, 1, double(bz)
close,1

end