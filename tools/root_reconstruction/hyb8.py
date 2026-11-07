import numpy as np,json,os
from scipy.ndimage import gaussian_filter
from PIL import Image
vol=np.load('hyb_vol2.npy').astype(np.float32);lab=np.load('hyb_lab3.npy');u1,u2,A,c=np.load('hyb_frame.npy')
R=np.load('hyb_Rfin.npy');HS=np.arange(-6,40.25,0.25);T=json.load(open('hyb_truth.json'))
TH=np.radians(np.arange(0,360,1.0))
rows=[int(np.argmin(abs(HS-h))) for h in np.arange(-2,1.01,0.25)]
r=R[rows].mean(0);F=np.fft.rfft(r);F[7:]=0;r=np.fft.irfft(F,len(r))
x=r*np.cos(TH);y=r*np.sin(TH);xs,ys=np.roll(x,-1),np.roll(y,-1);cr=x*ys-xs*y;Ar=cr.sum()/2
cx=np.sum((x+xs)*cr)/(6*Ar);cy=np.sum((y+ys)*cr)/(6*Ar)
X_,Y_=x-cx,y-cy;Xs,Ys=np.roll(X_,-1),np.roll(Y_,-1);a=X_*Ys-Xs*Y_
Iyy=np.sum(a*(Y_**2+Y_*Ys+Ys**2))/12/Ar;Ixx=np.sum(a*(X_**2+X_*Xs+Xs**2))/12/Ar;Ixy=np.sum(a*(X_*Ys+2*X_*Y_+2*Xs*Ys+Xs*Y_))/24/Ar
w,V=np.linalg.eigh(np.array([[Ixx,Ixy],[Ixy,Iyy]]));sa,sb=2*np.sqrt(w[1]),2*np.sqrt(w[0]);ea=V[:,1]
print('fit semi-axes',sa,sb,'centroid',cx,cy,'contour area',Ar)
# keep area of the real contour: scale ellipse so area matches
k=np.sqrt(Ar/(np.pi*sa*sb));sa*=k;sb*=k
eA=ea[0]*u1+ea[1]*u2;eB=np.cross(A,eA);eB/=np.linalg.norm(eB)
cc=c+cx*u1+cy*u2
N=220;VS=.5;ax=((np.arange(N)+.5)*VS).astype(np.float32);Z,Y,X=np.meshgrid(ax,ax,ax,indexing='ij')
dX,dY,dZ=X-cc[0],Y-cc[1],Z-cc[2];H=dX*A[0]+dY*A[1]+dZ*A[2];P=dX*eA[0]+dY*eA[1]+dZ*eA[2];Q=dX*eB[0]+dY*eB[1]+dZ*eB[2]
del dX,dY,dZ,X,Y,Z
# LVOT tapers slightly below the annulus
scale=np.where(H<0,1-0.03*np.clip(-H/6,0,1),1.0)
rho=np.sqrt((P/(sa*scale))**2+(Q/(sb*scale))**2)
dist=(rho-1)*np.minimum(sa,sb)*scale            # approx signed distance (mm)
wl=np.clip(np.minimum((H+7)/2.5,(2.5-H)/1.5),0,1)  # lumen raise band: h in [-7, 2.5]
wr=np.clip(np.minimum((H+6)/2.0,(1.5-H)/1.5),0,1)  # rim band: h in [-6, 1.5]
lum=T['lumenHU'];prot=(lab==3)|(lab==10)|(lab==8)|(lab==9)
inside=np.clip(0.5-dist/0.6,0,1)
rim=np.clip(0.5-np.abs(dist-1.4)/0.9,0,1)
sm=gaussian_filter(vol,2.0)
v2=vol+inside*wl*np.clip(lum-10-sm,0,None)
v2=np.where(prot,vol,v2)
sm2=gaussian_filter(v2,2.0)
v2=v2-0.8*rim*wr*np.clip(sm2-130,0,None)
v2=np.where(prot,vol,v2)
v2=np.clip(v2,-1024,3071).astype(np.int16)
lab[(dist<0)&(np.abs(H)<0.75)&~prot]=5
np.save('hyb_vol3.npy',v2);np.save('hyb_lab4.npy',lab)
per=float(np.pi*(3*(sa+sb)-np.sqrt((3*sa+sb)*(sa+3*sb))))
T.update(dict(max=round(2*sa,1),min=round(2*sb,1),area=round(float(np.pi*sa*sb)),per=round(per,1),c=[round(float(q),2) for q in cc]))
json.dump(T,open('hyb_truth.json','w'));print({k:T[k] for k in('max','min','area','per','c','n_toLV')})
u=(v2.astype(np.int32)+1024).astype(np.uint16);files=[]
for i,z0 in enumerate(range(0,N,64)):
    z1=min(z0+64,N);c_=u[z0:z1].reshape(-1,N);ll=lab[z0:z1].reshape(-1,N)
    rgb=np.stack([(c_>>8).astype(np.uint8),(c_&255).astype(np.uint8),ll],-1);fn=f'/home/claude/ctp/data/root_{i:02d}.png'
    Image.fromarray(rgb,'RGB').save(fn,optimize=True,compress_level=9);files.append([z0,z1,os.path.getsize(fn)])
print(files)
