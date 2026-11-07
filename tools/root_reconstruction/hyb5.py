import numpy as np,json
from scipy.ndimage import map_coordinates,median_filter,gaussian_filter1d
up=np.load('hyb_up.npy').astype(np.float32);lab=np.load('hyb_lab.npy');u1,u2,A,c=np.load('hyb_frame.npy')
def sampv(vol,P):return map_coordinates(vol,[P[...,2]/.5-.5,P[...,1]/.5-.5,P[...,0]/.5-.5],order=1,cval=-1000)
TH=np.radians(np.arange(0,360,2.0));HS=np.arange(-14,44.5,0.5);RS=np.arange(0,30,0.25)
d=np.cos(TH)[:,None,None]*u1+np.sin(TH)[:,None,None]*u2
Rimg=np.zeros((len(HS),len(TH)))
for i,h in enumerate(HS):
    P=c+A*h+d*RS[None,:,None];val=sampv(up,P);b=val<210;b[:,0]=False;Rimg[i]=np.where(b.any(1),RS[b.argmax(1)],RS[-1])
# annulus label outline in plane coordinates
z,y,x=np.nonzero(lab==5);P=(np.stack([x,y,z],1)+.5)*.5;q=P-c
hh=q@A;a1=q@u1;a2=q@u2;th=np.arctan2(a2,a1);rr=np.hypot(a1,a2)
bins=((np.degrees(th)%360)//2).astype(int);Rann=np.zeros(180)
for k in range(180):
    s=rr[bins==k];Rann[k]=np.percentile(s,90) if len(s)>3 else np.nan
ok=~np.isnan(Rann);Rann=np.interp(np.arange(180),np.arange(180)[ok],Rann[ok],period=180)
Rann=gaussian_filter1d(median_filter(Rann,7,mode='wrap'),2,mode='wrap')
print('annulus label: plane offsets h mean/sd',hh.mean().round(2),hh.std().round(2),'R mean',Rann.mean().round(1),'min',Rann.min().round(1),'max',Rann.max().round(1))
# capped sinus radius: cap leakage per height at 1.12 x the 50th percentile, never below annulus-based floor
R=np.zeros_like(Rimg)
for i,h in enumerate(HS):
    r=Rimg[i];cap=np.percentile(r,50)*1.12;R[i]=np.minimum(r,cap)
R=median_filter(R,(3,9),mode='wrap');R=gaussian_filter1d(gaussian_filter1d(R,2,axis=1,mode='wrap'),2,axis=0)
for h in (-3,0,3,6,10,14,18,22,26,30,36):
    i=np.argmin(abs(HS-h));print(f'h={h:4} capped meanR={R[i].mean():5.1f} min {R[i].min():5.1f} max {R[i].max():5.1f}')
np.save('hyb_Rcap.npy',R);np.save('hyb_Rann.npy',Rann)
# annulus dims from label outline
x_=Rann*np.cos(TH);y_=Rann*np.sin(TH);xs,ys=np.roll(x_,-1),np.roll(y_,-1);cr=x_*ys-xs*y_
area=abs(cr.sum())/2;per=np.sum(np.hypot(xs-x_,ys-y_))
cx=np.sum((x_+xs)*cr)/(3*cr.sum());cy=np.sum((y_+ys)*cr)/(3*cr.sum())
X_,Y_=x_-cx,y_-cy;Xs,Ys=np.roll(X_,-1),np.roll(Y_,-1);a=X_*Ys-Xs*Y_;Ar=a.sum()/2
Ixx=np.sum(a*(Y_**2+Y_*Ys+Ys**2))/12;Iyy=np.sum(a*(X_**2+X_*Xs+Xs**2))/12;Ixy=np.sum(a*(X_*Ys+2*X_*Y_+2*Xs*Ys+Xs*Y_))/24
w,V=np.linalg.eigh(np.array([[Iyy,Ixy],[Ixy,Ixx]])/Ar)
print('annulus area',round(area),'perimeter',round(per,1),'Dmax',round(4*np.sqrt(w[1]),1),'Dmin',round(4*np.sqrt(w[0]),1),'centroid',round(cx,2),round(cy,2))
json.dump(dict(area=float(area),per=float(per),dmax=float(4*np.sqrt(w[1])),dmin=float(4*np.sqrt(w[0])),cx=float(cx),cy=float(cy)),open('hyb_ann.json','w'))
