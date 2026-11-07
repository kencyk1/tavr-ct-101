import numpy as np,json
from scipy.ndimage import map_coordinates,median_filter,gaussian_filter1d,gaussian_filter
up=np.load('hyb_up.npy').astype(np.float32);sh=np.load('hyb_sh.npy').astype(np.float32);lab=np.load('hyb_lab.npy')
u1,u2,A,c=np.load('hyb_frame.npy');N=220;VS=.5
def sampv(vol,P):return map_coordinates(vol,[P[...,2]/.5-.5,P[...,1]/.5-.5,P[...,0]/.5-.5],order=1,cval=-1000)
DT=1.0;TH=np.radians(np.arange(0,360,DT));HS=np.arange(-6,40.25,0.25);RS=np.arange(0,30,0.2)
d=np.cos(TH)[:,None,None]*u1+np.sin(TH)[:,None,None]*u2
Rimg=np.zeros((len(HS),len(TH)))
for i,h in enumerate(HS):
    P=c+A*h+d*RS[None,:,None];val=sampv(up,P);b=val<250;b[:,0]=False;r=np.where(b.any(1),RS[b.argmax(1)],RS[-1])
    Rimg[i]=np.minimum(r,np.percentile(r,50)*1.12)
R=median_filter(Rimg,(5,15),mode='wrap');R=gaussian_filter1d(gaussian_filter1d(R,4,axis=1,mode='wrap'),3,axis=0)
iH=lambda h:int(np.argmin(abs(HS-h)))
# sinus lobes at h=10..14
prof=R[iH(10):iH(14)].mean(0);pk=[]
for k in range(len(TH)):
    if prof[k]==max(prof[(k+j)%len(TH)] for j in range(-25,26)):pk.append(k)
print('lobe peaks deg',[int(np.degrees(TH[k])) for k in pk],'R',[round(prof[k],1) for k in pk])
np.save('hyb_Rfin.npy',R)
# cusp centres from lobes and commissures midway
cen={'RCC':-23.0,'LCC':118.0,'NCC':-105.0}
order=['RCC','LCC','NCC']
def mid(a,b):
    d=((b-a)%360);return a+d/2
comm={('RCC','LCC'):mid(-23,118),('LCC','NCC'):mid(118,255),('NCC','RCC'):mid(255,337)}
print('commissures',comm)
Hc,Hf,T,HU_LEAF=20.0,16.0,2.2,95.0
ax=((np.arange(N)+.5)*VS).astype(np.float32)
Z,Y,X=np.meshgrid(ax,ax,ax,indexing='ij')
dX,dY,dZ=X-c[0],Y-c[1],Z-c[2]
H=dX*A[0]+dY*A[1]+dZ*A[2];Q1=dX*u1[0]+dY*u1[1]+dZ*u1[2];Q2=dX*u2[0]+dY*u2[1]+dZ*u2[2]
Rr=np.hypot(Q1,Q2);THd=np.degrees(np.arctan2(Q2,Q1))%360
del dX,dY,dZ,Q1,Q2
def Rlook(h,thd):
    return map_coordinates(R,[np.clip((h-HS[0])/0.25,0,len(HS)-1),(thd%360)/DT],order=1,mode='grid-wrap')
roi=(H>-1.5)&(H<Hc+1)&(Rr<26)
idx=np.nonzero(roi);h=H[idx];th=THd[idx];r=Rr[idx]
sd=np.full(h.shape,99.0,np.float32);sector=np.full(h.shape,-1)
bounds={'RCC':(comm[('NCC','RCC')]-360,comm[('RCC','LCC')]),'LCC':(comm[('RCC','LCC')],comm[('LCC','NCC')]),'NCC':(comm[('LCC','NCC')]-360,comm[('NCC','RCC')]-360)}
leaf_par={}
for i,nm in enumerate(order):
    lo,hi=bounds[nm];ctr=cen[nm]
    dd=((th-ctr+180)%360)-180
    inside=(dd>=lo-ctr)&(dd<=hi-ctr)
    half=np.where(dd<0,ctr-lo,hi-ctr);t=np.clip(np.abs(dd)/half,0,1)
    hatt=Hc*(1-np.cos(np.pi/2*t));ok=inside&(hatt<Hf-1.0)
    s=(h-hatt)/np.maximum(Hf-hatt,1.0);ok&=(s>=0)&(s<=1)
    ratt=Rlook(hatt,th);rfree=4.5+1.0*np.cos(np.pi*t)
    sc=np.clip(s,0,1)**0.75;rl=ratt*(1-sc)+rfree*sc
    sd=np.where(ok,np.minimum(sd,np.abs(r-rl)-T/2),sd)
    leaf_par[nm]=(lo,hi,ctr)
a=np.clip(0.5-sd/0.4,0,1).astype(np.float32)
aL=np.zeros(sh.shape,np.float32);aL[idx]=a;aL=np.clip(gaussian_filter(aL,0.85)*1.15,0,1)
vol=sh*(1-aL)+HU_LEAF*aL
lab2=lab.copy();lab2[aL>0.35]=3
aC=np.zeros(sh.shape,np.float32);hC=np.zeros(sh.shape,np.float32)
# leaflet calcium nodules
rng=np.random.default_rng(5)
def leaf_pt(nm,frac_ang,s):
    lo,hi,ctr=leaf_par[nm];half=(ctr-lo) if frac_ang<0 else (hi-ctr);thd=ctr+frac_ang*half;t=abs(frac_ang)
    hatt=Hc*(1-np.cos(np.pi/2*t));hh=hatt+s*max(Hf-hatt,1);ra=float(Rlook(np.array([hatt]),np.array([thd]))[0]);rf=4.5+np.cos(np.pi*t)
    rr_=ra*(1-s**.75)+rf*s**.75;tt=np.radians(thd);return c+A*hh+rr_*(np.cos(tt)*u1+np.sin(tt)*u2)
nod=[]
for nm,n in (('NCC',3),('LCC',2),('RCC',2)):
    for _ in range(n):
        p=leaf_pt(nm,rng.uniform(-.6,.6),rng.uniform(.15,.7));r0=rng.uniform(1.0,1.8);nod.append((p,r0))
        for _ in range(2):nod.append((p+rng.standard_normal(3)/np.sqrt(3)*r0*.8,r0*rng.uniform(.4,.7)))
for p,r0 in nod:
    lo_=np.maximum(((p-r0-2)/VS).astype(int),0);hi_=np.minimum(((p+r0+2)/VS).astype(int)+1,N)
    sl=(slice(lo_[2],hi_[2]),slice(lo_[1],hi_[1]),slice(lo_[0],hi_[0]))
    k=np.sqrt((X[sl]-p[0])**2+(Y[sl]-p[1])**2+(Z[sl]-p[2])**2)/r0
    al=np.clip(0.5-(k-1)*r0/0.35,0,1);hu=rng.uniform(900,1200);aC[sl]=np.maximum(aC[sl],al);hC[sl]=np.where(al>0,hu,hC[sl])
aC=np.clip(gaussian_filter(aC,0.9)*1.2,0,1);hC=np.maximum(gaussian_filter(hC,0.9)/np.maximum(gaussian_filter((hC>0).astype(np.float32),0.9),1e-3),0)
vol=np.maximum(vol,vol*(1-aC)+hC*aC);lab2[aC>0.4]=10
vol=np.clip(vol,-1024,3071).astype(np.int16)
np.save('hyb_vol.npy',vol);np.save('hyb_lab2.npy',lab2)
# annulus truth from image contour at h=0
r0=R[iH(0)];F=np.fft.rfft(r0);F[7:]=0;r0=np.fft.irfft(F,len(r0));x_=r0*np.cos(TH);y_=r0*np.sin(TH);xs,ys=np.roll(x_,-1),np.roll(y_,-1);cr=x_*ys-xs*y_
area=abs(cr.sum())/2;per=float(np.sum(np.hypot(xs-x_,ys-y_)));cx=np.sum((x_+xs)*cr)/(3*cr.sum());cy=np.sum((y_+ys)*cr)/(3*cr.sum())
X_,Y_=x_-cx,y_-cy;Xs,Ys=np.roll(X_,-1),np.roll(Y_,-1);aa=X_*Ys-Xs*Y_;Ar=aa.sum()/2
Ixx=np.sum(aa*(Y_**2+Y_*Ys+Ys**2))/12;Iyy=np.sum(aa*(X_**2+X_*Xs+Xs**2))/12;Ixy=np.sum(aa*(X_*Ys+2*X_*Y_+2*Xs*Ys+Xs*Y_))/24
w,V=np.linalg.eigh(np.array([[Iyy,Ixy],[Ixy,Ixx]])/Ar)
ref_c=c+cx*u1+cy*u2
truth=dict(area=round(float(area)),per=round(per,1),max=round(float(4*np.sqrt(w[1])),1),min=round(float(4*np.sqrt(w[0])),1),c=[round(float(q),2) for q in ref_c],n_toLV=[round(float(-q),4) for q in A],
 sinus_d=round(float(2*R[iH(12)].mean()),1),stj_d=round(float(2*R[iH(26)].mean()),1))
json.dump(truth,open('hyb_truth.json','w'));print(truth)
