pa=dblarr(2001)
pa_t=dblarr(2001)
mu=dblarr(2001)
openr,1,'lowlou_y.txt'
readf,1,pa
close,1

openr,1,'lowlou_dy.txt'
readf,1,pa_t
close,1

openr,1,'lowlou_x.txt'
readf,1,mu
close,1


l1=-0.525
l2=0.525

dd=0.015
n=1
; --------------------------------------------------------------------
;The eq is solved by fix P' not coefficent \alpha, 
;in order to have the same solution solution with P'(mu=-1;mu=1)=10,
;we change it as follows.
alpha=1*(pa_t[0]/10.0)^(1./n)
pa=pa*(10./pa_t[0])
pa_t=pa_t*(10/pa_t[0])
; --------------------------------------------------------------------
dn=fix((l2-l1)/dd)
x0=fltarr(dn,dn,dn)
y0=fltarr(dn,dn,dn)
z0=fltarr(dn,dn,dn)

for i=0,dn-1 do x0[i,*,*]=i*dd
for i=0,dn-1 do y0[*,i,*]=l1+i*dd
for i=0,dn-1 do z0[*,*,i]=l1+i*dd

rr=sqrt(x0^2+y0^2+z0^2)
mu0=z0/rr
pa_t1=Interpol(pa_t,mu,mu0)
pa_1=Interpol(pa,mu,mu0)

cost=z0/rr
sint=sqrt(1-cost^2)
sinp=y0/sqrt(x0^2+y0^2)
cosp=x0/sqrt(x0^2+y0^2)

br=-1*pa_t1/(rr^(2+n))
bt=(n/(rr^(2+n)*sint))*pa_1
bp = alpha*(1/(rr^(2+n)*sint))*pa_1*abs(pa_1)^(1.0/n)

bt(where(sint eq 0 or rr eq 0))=0.0
bp(where(sint eq 0 or rr eq 0))=0.0

bx = br*sint*cosp + bt*cost*cosp - bp*sinp
by = br*sint*sinp + bt*cost*sinp + bp*cosp
bz = br*cost          - (n/(rr^(2+n)))*pa_1

bb=sqrt(bx^2+by^2+bz^2)
bbr=sqrt(br^2+bt^2+bp^2)

jx = 0.5*(shift(by,[0,0,-1])-shift(by,[0,0,1])) - 0.5*(shift(bz,[0,-1,0])-shift(bz,[0,1,0])) 
jy = 0.5*(shift(bx,[0,0,-1])-shift(bx,[0,0,1])) - 0.5*(shift(bz,[-1,0,0])-shift(bz,[1,0,0])) 
jz = 0.5*(shift(by,[-1,0,0])-shift(by,[1,0,0])) - 0.5*(shift(bx,[0,-1,0])-shift(bx,[0,1,0])) 

mark=fix(1.0*(rr gt 0.2))
bx=bx*mark
by=by*mark
bz=bz*mark
siz=size(bz)
alpha=siz[1]*(jx*bx+jy*by+jz*bz)/(bx^2+by^2+bz^2)

; bx0=reform(bx[1:128,40,1:128])
; by0=reform(by[1:128,40,1:128])
; bz0=reform(bz[1:128,40,1:128])
; alpha0=reform(alpha[1:128,40,1:128])

;bx0=reform(bx[40,1:128,1:128])
;by0=reform(by[40,1:128,1:128])
;bz0=reform(bz[40,1:128,1:128])
;alpha0=reform(alpha[40,1:128,1:128])


;bx1=rebin(bx0,64,64)
;by1=rebin(by0,64,64)
;bz1=rebin(bz0,64,64)
;alpha1=rebin(alpha0,64,64)


bx1=bx[20,1:64,1:64]/max(abs(bx[20,1:64,1:64]))
by1=by[20,1:64,1:64]/max(abs(bx[20,1:64,1:64]))
bz1=bz[20,1:64,1:64]/max(abs(bx[20,1:64,1:64]))
alpha1=alpha[20,1:64,1:64]

OPENW, 1, './data/bz0.dat'
WRITEU, 1, double(bx1)
close,1

OPENW, 1, './data/alpha0.dat'
WRITEU, 1,double(alpha1)
close,1

OPENW, 1, './data/sig0.dat'
WRITEU, 1, double(abs(bx1))
close,1



OPENW, 1, './bz_ll.dat'
WRITEU, 1, double(bx[20:69,1:64,1:64]/max(abs(bx[20,1:64,1:64])))
close,1

OPENW, 1, './bx_ll.dat'
WRITEU, 1,double(by[20:69,1:64,1:64]/max(abs(bx[20,1:64,1:64])))
close,1

OPENW, 1, './by_ll.dat'
WRITEU, 1, double(bz[20:69,1:64,1:64]/max(abs(bx[20,1:64,1:64])))
close,1






; plot_image,alog10(reform(frb[150,*,*])),max=5,min=0
; plot_image,reform(frb[150,*,*]),max=100,min=-100

; vec_2vtk,bx,by,bz,'LowLou_n1_a1.vtk'

; bxyz=fltarr(200,200,200,3)
; jxyz=fltarr(200,200,200,3)
; bxyz[*,*,*,0]=float(bx)
; bxyz[*,*,*,1]=float(by)
; bxyz[*,*,*,2]=float(bz)

; jxyz[*,*,*,0]=float(jx)
; jxyz[*,*,*,1]=float(jy)
; jxyz[*,*,*,2]=float(jz)
; data={bxyz:bxyz,jxyz:jxyz,dim:3}
; writevtk,'lowlou_n1_a1.vtk',data
end