import numpy as np,json
from scipy.ndimage import map_coordinates,gaussian_filter,zoom,gaussian_filter1d
from skimage.graph import MCP_Geometric
up=np.load('hyb_up.npy').astype(np.float32);u1,u2,A,c=np.load('hyb_frame.npy')
R=np.load('hyb_Rfin.npy');HS=np.arange(-6,40.25,0.25)
def Rw(h,deg):return float(map_coordinates(R,[[np.clip((h-HS[0])/0.25,0,len(HS)-1)],[deg%360]],order=1,mode='grid-wrap')[0])
v1=gaussian_filter(up,0.8)[::2,::2,::2]          # 1 mm grid, index (z,y,x), centre of voxel k at (2k+1)*0.5 mm
n=v1.shape[0]
cost=np.full(v1.shape,6.0,np.float32)
cost[(v1>-220)&(v1<-25)]=1.0
cost[v1>=150]=80.0; cost[v1>=600]=4.0; cost[v1<-400]=60.0
# forbid the root lumen interior
ax=(np.arange(n)+.5)*1.0
Z,Y,X=np.meshgrid(ax,ax,ax,indexing='ij');d=np.stack([X-c[0],Y-c[1],Z-c[2]],-1)
H=d@A;th=np.degrees(np.arctan2(d@u2,d@u1));rr=np.hypot(d@u1,d@u2)
Rm=map_coordinates(R,[np.clip((H-HS[0])/0.25,0,len(HS)-1),th%360],order=1,mode='grid-wrap')
inroot=(H>-6)&(H<40)&(rr<Rm+1.0);cost[inroot]=500.0
def idx(p):return tuple(np.clip(np.round(np.asarray(p)[::-1]/1.0-.5).astype(int),0,n-1))
def wpt(i):return (np.array(i[::-1])+.5)*1.0
def route(start,dirv,dmin,dmax,name,CUMK=2.2,extra=None):
    m=MCP_Geometric(cost);cum,_=m.find_costs([idx(start)])
    E=np.sqrt((X-start[0])**2+(Y-start[1])**2+(Z-start[2])**2)
    proj=((X-start[0])*dirv[0]+(Y-start[1])*dirv[1]+(Z-start[2])*dirv[2])
    cand=(E>dmin)&(E<dmax)&(cost<=1.0)&(cum<CUMK*E)&(X>4)&(X<n-4)&(Y>4)&(Y<n-4)&(Z>4)&(Z<n-4)
    if extra is not None: cand&=extra(X,Y,Z,start)
    if not cand.any(): print(name,'no candidate');return None
    score=np.where(cand,proj/E-0.004*cum,-1e9);k=np.unravel_index(np.argmax(score),score.shape)
    path=np.array([wpt(i) for i in m.traceback(k)],float)
    path=np.vstack([path[:1],path]);sm=np.stack([gaussian_filter1d(path[:,j],2.5,mode='nearest') for j in range(3)],1)
    L=np.sum(np.linalg.norm(np.diff(sm,axis=0),axis=1));hu=map_coordinates(v1,[sm[:,2]-.5,sm[:,1]-.5,sm[:,0]-.5],order=1)
    print(name,'len',round(L,1),'end',k,'cum/E',round(float(cum[k]/E[k]),2),'HU along median',round(float(np.median(hu)),0),'frac>150',round(float((hu>150).mean()),2))
    return sm
def ost(deg,h,extra):
    t=np.radians(deg);rad=np.cos(t)*u1+np.sin(t)*u2;rw=Rw(h,deg);return c+A*h+rad*(rw+extra),c+A*h+rad*(rw-1.0),rad
lm_s,lm_in,lm_rad=ost(106,12.1,2.5);rca_s,rca_in,rca_rad=ost(-12,20.3,2.5)
nv=lambda v:np.asarray(v,float)/np.linalg.norm(v)
LMp=route(lm_s,nv([0.8,-0.2,0.2]),8,11,'LM')
B=LMp[-1]
LADp=route(B,nv([0.55,-0.7,0.45]),14,24,'LAD',CUMK=2.6,extra=lambda X,Y,Z,s:((Y-s[1])<-6)&((X-s[0])>0))
LCXp=route(B,nv([0.45,0.75,0.45]),14,26,'LCx',CUMK=2.6,extra=lambda X,Y,Z,s:((Y-s[1])>4))
LAD=np.vstack([LMp,LADp[1:]]);LCX=np.vstack([LMp,LCXp[1:]])
np.save('lm_bif.npy',B);print('LM length',round(float(np.sum(np.linalg.norm(np.diff(LMp,axis=0),axis=1))),1))
seed,_,_=ost(-52,19.0,3.5)
print('rca seed HU',float(map_coordinates(v1,[[seed[2]-.5],[seed[1]-.5],[seed[0]-.5]],order=1)[0]),'cost',cost[idx(seed)])
RCAr=route(seed,nv([-0.55,-0.45,0.7]),18,30,'RCA',CUMK=2.6,extra=lambda X,Y,Z,s:((X-s[0])<-2))
arc=[];o_in=None
for k,deg in enumerate(np.linspace(-35,-52,8)):
    p,pin,_=ost(deg,19.0,2.0+1.5*k/7);arc.append(p)
    if k==0:o_in=pin
RCA=np.vstack([np.array(arc),RCAr[1:]]);RCA=np.stack([gaussian_filter1d(RCA[:,j],1.5,mode='nearest') for j in range(3)],1)
rca_in=o_in;rca_s=arc[0]
np.savez('routes.npz',LAD=LAD,LCX=LCX,RCA=RCA,lm_in=lm_in,rca_in=rca_in,lm_s=lm_s,rca_s=rca_s)
