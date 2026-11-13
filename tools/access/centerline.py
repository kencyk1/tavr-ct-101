import numpy as np, json
from scipy.ndimage import label as cc, binary_dilation, gaussian_filter1d, map_coordinates, gaussian_filter
vol=np.load('vol.npy');lab=np.load('lab.npy')
VS=1.5
def comps(mask):
    l,n=cc(mask);return [(l==k) for k in range(1,n+1)]
def cen(m): yy,xx=np.nonzero(m);return np.array([xx.mean(),yy.mean()])
def track(L,z0,z1,zend):
    pts=[]
    # aorta above the bifurcation
    for z in range(z0,z1):
        m=(lab[z]==1)&(vol[z]>180)
        if not m.any(): m=lab[z]==1
        pts.append((z,*cen(m)))
    prev=np.array(pts[-1][1:]);vel=np.zeros(2)
    for z in range(z1,zend+1):
        cand=lab[z]==L
        lum=vol[z]>180
        if cand.any():
            best=None
            for c in comps(cand):
                cc_=cen(c);d=np.linalg.norm(cc_-(prev+vel))
                if best is None or d<best[0]:best=(d,c)
            m=best[1]
            m2=binary_dilation(m,iterations=2)&lum
            if m2.sum()>=4: m=m2
        else:
            yy,xx=np.mgrid[:lum.shape[0],:lum.shape[1]];p=prev+vel
            near=lum&((xx-p[0])**2+(yy-p[1])**2<9**2)
            if not near.any():break
            best=None
            for c in comps(near):
                cc_=cen(c);d=np.linalg.norm(cc_-p)
                if best is None or d<best[0]:best=(d,c)
            m=best[1]
        c_=cen(m);vel=0.6*vel+0.4*(c_-prev);prev=c_;pts.append((z,*c_))
    P=np.array(pts,float)  # z, x, y (voxels)
    W=np.stack([(P[:,1]+.5)*VS,(P[:,2]+.5)*VS,(P[:,0]+.5)*VS],1)  # world x,y,z (mm)
    # smooth and resample at 1 mm
    for k in range(3): W[:,k]=gaussian_filter1d(W[:,k],2.0,mode='nearest')
    s=np.concatenate([[0],np.cumsum(np.linalg.norm(np.diff(W,axis=0),axis=1))])
    S=np.arange(0,s[-1],1.0);R=np.stack([np.interp(S,s,W[:,k]) for k in range(3)],1)
    for k in range(3): R[:,k]=gaussian_filter1d(R[:,k],3.0,mode='nearest')
    return R
out={}
for side,L in (('R',7),('L',6)):
    R=track(L,322,352,496);out[side]=R.round(2).tolist();print(side,len(R),'mm long; end',R[-1].round(1))
json.dump(out,open('centerlines.json','w'))
