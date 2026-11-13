import numpy as np, json
from scipy.ndimage import map_coordinates, label as cc, binary_fill_holes
vol=np.load('vol.npy').astype(np.float32);C=json.load(open('centerlines.json'))
def frames(P):
    T=np.gradient(P,axis=0);T/=np.linalg.norm(T,axis=1,keepdims=True)
    n=np.cross(T[0],[0,1,0]);n/=np.linalg.norm(n);N=[n]
    for i in range(1,len(P)):
        v=N[-1]-T[i]*(N[-1]@T[i]);v/=np.linalg.norm(v);N.append(v)
    N=np.array(N);B=np.cross(T,N);return T,N,B
S=lambda Q:map_coordinates(vol,[Q[...,2]/1.5-.5,Q[...,1]/1.5-.5,Q[...,0]/1.5-.5],order=3,cval=-1000)
px=0.2;u=np.arange(-15,15,px)+px/2;UU,VV=np.meshgrid(u,u,indexing='xy')
angs=np.radians(np.arange(0,180,2))
res={}
for side in ('R','L'):
    P=np.array(C[side]);T,N,B=frames(P);rows=[]
    for i in range(len(P)):
        Q=P[i]+UU[...,None]*N[i]+VV[...,None]*B[i];im=S(Q)
        c=im[(np.abs(UU)<1.5)&(np.abs(VV)<1.5)];lum=float(np.median(c));thr=(lum+50)/2
        m=(im>thr)&(im<lum+300)
        l,n=cc(m);k=l[len(u)//2,len(u)//2]
        if k==0:
            d2=UU**2+VV**2;cand=[(d2[l==j].min(),j) for j in range(1,n+1)] if n else []
            k=min(cand)[1] if cand else 0
        if k==0: rows.append([i,0,0,0,lum,0]);continue
        mk=binary_fill_holes(l==k);yy,xx=UU[mk],VV[mk]
        area=mk.sum()*px*px;deq=2*np.sqrt(area/np.pi)
        w=[(np.ptp(yy*np.cos(a)+xx*np.sin(a))+px) for a in angs]
        calc=float(im[(UU**2+VV**2<(deq/2+2.5)**2)&~mk].max()) if True else 0
        rows.append([i,round(deq,2),round(min(w),2),round(max(w),2),round(lum),round(calc)])
    res[side]=rows
    A=np.array(rows)
    zz=P[:,2]/1.5-.5
    print(side)
    for i in range(0,len(P),10): print(f'  s{i:4d} z{zz[i]:6.1f}  deq {A[i,1]:5.1f} min {A[i,2]:5.1f} max {A[i,3]:5.1f} lum {A[i,4]:5.0f} wallmax {A[i,5]:5.0f}')
json.dump(res,open('diam.json','w'))
