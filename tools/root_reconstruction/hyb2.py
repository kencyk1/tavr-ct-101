import numpy as np
from scipy.ndimage import map_coordinates,median_filter,gaussian_filter1d
up=np.load('hyb_up.npy').astype(np.float32);u1,u2,A,c=np.load('hyb_frame.npy')
def samp(P):return map_coordinates(up,[P[...,2]/.5-.5,P[...,1]/.5-.5,P[...,0]/.5-.5],order=1,cval=-1000)
TH=np.radians(np.arange(0,360,2.0));HS=np.arange(-12,40.5,0.5);RS=np.arange(0,30,0.25)
thr=210
# recentre per height: iterate centroid of lumen ray ends
cent=np.zeros((len(HS),3))
R=np.zeros((len(HS),len(TH)))
for i,h in enumerate(HS):
    cc=c+A*h
    for it in range(3):
        d=np.cos(TH)[:,None,None]*u1+np.sin(TH)[:,None,None]*u2           # (T,1,3)
        P=cc+d*RS[None,:,None]
        val=samp(P)                                                       # (T,R)
        below=val<thr; below[:,0]=False
        idx=np.where(below.any(1),below.argmax(1),len(RS)-1)
        r=RS[idx];pts=cc+np.cos(TH)[:,None]*u1*r[:,None]+np.sin(TH)[:,None]*u2*r[:,None]
        # robust centroid in-plane
        rm=median_filter(r,7,mode='wrap')
        pts=cc+np.cos(TH)[:,None]*u1*rm[:,None]+np.sin(TH)[:,None]*u2*rm[:,None]
        cc=pts.mean(0)
    cent[i]=cc;R[i]=median_filter(r,9,mode='wrap')
np.save('hyb_R.npy',R);np.save('hyb_cent.npy',cent)
for i,h in enumerate(HS):
    if h%2==0: print(f'h={h:5.1f} meanR={R[i].mean():5.1f} minR={R[i].min():5.1f} maxR={R[i].max():5.1f} centershift={np.linalg.norm(cent[i]-(c+A*h)):.1f}')
