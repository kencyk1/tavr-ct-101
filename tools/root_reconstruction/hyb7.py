import numpy as np,json,os
from scipy.ndimage import map_coordinates,gaussian_filter
from PIL import Image
vol=np.load('hyb_vol.npy').astype(np.float32);lab=np.load('hyb_lab2.npy');u1,u2,A,c=np.load('hyb_frame.npy')
R=np.load('hyb_Rfin.npy');HS=np.arange(-6,40.25,0.25);T=json.load(open('hyb_truth.json'));c0=np.array(T['c'])
N=220;VS=.5;ax=((np.arange(N)+.5)*VS).astype(np.float32);Z,Y,X=np.meshgrid(ax,ax,ax,indexing='ij')
def Rw(h,deg):return float(map_coordinates(R,[[np.clip((h-HS[0])/0.25,0,len(HS)-1)],[deg%360]],order=1,mode='grid-wrap')[0])
COR={'LM':dict(deg=106,low=10.0,r=2.1,tilt=0.15,lab=8),'RCA':dict(deg=-12,low=18.5,r=1.8,tilt=-0.45,lab=9)}
aT=np.zeros(vol.shape,np.float32);out={}
for nm,k in COR.items():
    hc=k['low']+k['r'];t=np.radians(k['deg']);rad=np.cos(t)*u1+np.sin(t)*u2;rw=Rw(hc,k['deg'])
    p0=c0+A*hc+rad*(rw-2.5);p1=c0+A*hc+rad*(rw+4.0);p2=p1+(rad+A*k['tilt'])/np.linalg.norm(rad+A*k['tilt'])*8
    best=np.full(vol.shape,1e9,np.float32)
    for a_,b_ in ((p0,p1),(p1,p2)):
        v=b_-a_;wx,wy,wz=X-a_[0],Y-a_[1],Z-a_[2];tt=np.clip((wx*v[0]+wy*v[1]+wz*v[2])/v.dot(v),0,1)
        best=np.minimum(best,np.sqrt((wx-tt*v[0])**2+(wy-tt*v[1])**2+(wz-tt*v[2])**2))
    a=np.clip(0.5-(best-k['r'])/0.4,0,1);aT=np.maximum(aT,a);lab[(a>0.5)&(lab!=10)]=k['lab']
    out[nm]=dict(lower_edge_mm=k['low'],diameter_mm=2*k['r'],angle_deg=k['deg'])
aT=np.clip(gaussian_filter(aT,0.85)*1.1,0,1)
lum=np.percentile(vol[(lab==2)],60)
vol=np.where(vol>lum,vol,vol*(1-aT)+lum*aT)
vol=np.clip(vol,-1024,3071).astype(np.int16)
np.save('hyb_vol2.npy',vol);np.save('hyb_lab3.npy',lab)
T['coronaries']=out;T['lumenHU']=float(lum);json.dump(T,open('hyb_truth.json','w'));print(T)
u=(vol.astype(np.int32)+1024).astype(np.uint16);files=[]
for i,z0 in enumerate(range(0,N,64)):
    z1=min(z0+64,N);cc=u[z0:z1].reshape(-1,N);ll=lab[z0:z1].reshape(-1,N)
    rgb=np.stack([(cc>>8).astype(np.uint8),(cc&255).astype(np.uint8),ll],-1);fn=f'/home/claude/ctp/data/root_{i:02d}.png'
    Image.fromarray(rgb,'RGB').save(fn,optimize=True,compress_level=9);files.append([z0,z1,os.path.getsize(fn)])
print(files,sum(f[2] for f in files)/1e6,'MB');print(np.unique(lab))
