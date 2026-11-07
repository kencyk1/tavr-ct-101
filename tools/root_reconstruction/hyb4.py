import numpy as np,json
from scipy.ndimage import map_coordinates,median_filter,gaussian_filter,gaussian_filter1d
up=np.load('hyb_up.npy').astype(np.float32);sh=np.load('hyb_sh.npy').astype(np.float32)
m=np.load('hyb_meta.npy');O=m[:3];u1,u2,A,c=np.load('hyb_frame.npy')
lb=np.load('lb.npy');N=220;VS=.5
ax=((np.arange(N)+.5)*VS).astype(np.float32)
Z,Y,X=np.meshgrid(ax+O[2],ax+O[1],ax+O[0],indexing='ij')
lab=map_coordinates(lb,[Z/1.5-.5,Y/1.5-.5,X/1.5-.5],order=0,cval=0).astype(np.uint8)
del Z,Y,X
def sampv(vol,P,order=1):return map_coordinates(vol,[P[...,2]/.5-.5,P[...,1]/.5-.5,P[...,0]/.5-.5],order=order,cval=-1000)
TH=np.radians(np.arange(0,360,2.0));HS=np.arange(-14,44.5,0.5);RS=np.arange(0,30,0.25)
d=np.cos(TH)[:,None,None]*u1+np.sin(TH)[:,None,None]*u2
inlab=np.isin(lab,[1,2,3,5]).astype(np.float32)
Rimg=np.zeros((len(HS),len(TH)));Rlab=np.zeros_like(Rimg)
for i,h in enumerate(HS):
    P=c+A*h+d*RS[None,:,None]
    val=sampv(up,P);b=val<210;b[:,0]=False;Rimg[i]=np.where(b.any(1),RS[b.argmax(1)],RS[-1])
    lv=sampv(inlab,P,0)<.5;lv[:,0]=False;Rlab[i]=np.where(lv.any(1),RS[lv.argmax(1)],0)
R=np.where(Rlab>0,np.minimum(Rimg,Rlab+1.0),Rimg)
R=median_filter(R,(3,9),mode='wrap');R=gaussian_filter1d(gaussian_filter1d(R,2,axis=1,mode='wrap'),1.5,axis=0)
np.save('hyb_Rc.npy',R)
i0=np.argmin(abs(HS-0))
for h in (-6,-3,0,3,6,10,14,18,22,26,30,36):
    i=np.argmin(abs(HS-h));print(f'h={h:4} meanR={R[i].mean():5.1f} min {R[i].min():5.1f} max {R[i].max():5.1f}')
# annulus truth from contour at h=0
r0=R[i0];pts=np.stack([r0*np.cos(TH),r0*np.sin(TH)],1)
x,y=pts[:,0],pts[:,1];xs,ys=np.roll(x,-1),np.roll(y,-1)
area=0.5*abs(np.sum(x*ys-xs*y));per=np.sum(np.hypot(xs-x,ys-y))
cx=np.sum((x+xs)*(x*ys-xs*y))/(6*area*np.sign(np.sum(x*ys-xs*y)));cy=np.sum((y+ys)*(x*ys-xs*y))/(6*area*np.sign(np.sum(x*ys-xs*y)))
# second moments of the polygon (filled) via triangle fan
def moments(x,y):
    a=x*np.roll(y,-1)-np.roll(x,-1)*y;A_=a.sum()/2
    Ixx=np.sum(a*(y**2+y*np.roll(y,-1)+np.roll(y,-1)**2))/12;Iyy=np.sum(a*(x**2+x*np.roll(x,-1)+np.roll(x,-1)**2))/12
    Ixy=np.sum(a*(x*np.roll(y,-1)+2*x*y+2*np.roll(x,-1)*np.roll(y,-1)+np.roll(x,-1)*y))/24
    return A_,Ixx,Iyy,Ixy
A_,Ixx,Iyy,Ixy=moments(x-cx,y-cy);C=np.array([[Iyy,Ixy],[Ixy,Ixx]])/A_;w,V=np.linalg.eigh(C)
dmax,dmin=4*np.sqrt(w[1]),4*np.sqrt(w[0])
print('annulus area',round(area),'per',round(per,1),'ellipse D',round(dmax,1),round(dmin,1),'centroid offset',round(cx,1),round(cy,1))
np.save('hyb_lab.npy',lab)
json.dump(dict(area=float(area),per=float(per),dmax=float(dmax),dmin=float(dmin),cx=float(cx),cy=float(cy)),open('hyb_ann.json','w'))
