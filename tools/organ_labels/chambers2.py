import numpy as np
from scipy.ndimage import gaussian_filter,binary_dilation,binary_erosion,binary_opening,distance_transform_edt,label as cc
from skimage.segmentation import watershed
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
v=np.load(S+'/v.npy').astype(np.float32);lab=np.load('labels_merged.npy')
heart=(lab==22)|(lab==4)
zz,yy,xx=np.where(heart);z0,z1,y0,y1,x0,x1=zz.min()-4,zz.max()+5,yy.min()-4,yy.max()+5,xx.min()-4,xx.max()+5
sl=(slice(z0,z1),slice(y0,y1),slice(x0,x1));O=np.array([z0,y0,x0])
V=gaussian_filter(v[sl],0.8);L=lab[sl];Hm=heart[sl]
cav=Hm&(V>=145)
near=lambda ids,it:binary_dilation(np.isin(L,ids),iterations=it)
mk=np.zeros(L.shape,np.int32)
mk[binary_erosion(L==4,iterations=2)&cav]=1
mk[near([48,49],3)&cav&~(L==4)]=2
mk[near([40,41,66],3)&cav&~(L==4)]=3
ws=watershed(-V,markers=mk,mask=cav,compactness=0.001)
RH=ws==3
def cen(m,which):
    z=np.where(m.any(axis=(1,2)))[0];zs=z[-4:] if which=='low' else z[:4];pts=np.argwhere(m&np.isin(np.arange(m.shape[0])[:,None,None],zs));return pts.mean(0)
svc=cen(L==40,'low');ivc=cen((L==41)&(np.arange(L.shape[0])[:,None,None]<=np.where(RH.any(axis=(1,2)))[0].max()+6),'high')
pa=cen(L==66,'low')
P=np.argwhere(RH);score=-P[:,1]*1.0+P[:,2]*0.6+P[:,0]*0.8   # anterior, leftward, inferior
apex=P[np.argmax(score)].astype(float)
print('svc',svc.round(1),'ivc',ivc.round(1),'pa',pa.round(1),'apex',apex)
def dseg(p,a,b):
    ab=b-a;t=np.clip(((p-a)@ab)/(ab@ab),0,1);return np.linalg.norm(p-(a+t[:,None]*ab),axis=1)
dRA=dseg(P.astype(float),svc,ivc);dRV=dseg(P.astype(float),pa,apex)
isRV=dRV<dRA*0.95
out=np.zeros(L.shape,np.uint8);out[ws==1]=4;out[ws==2]=96
pR=P[isRV];pA=P[~isRV];out[tuple(pR.T)]=95;out[tuple(pA.T)]=94
# smooth RA/RV boundary
for i in (94,95):
    m=binary_opening(out==i,iterations=1);out[(out==i)&~m]=0
# fill holes left by opening using nearest of RA/RV
hole=RH&(out==0)
if hole.any():
    d94=distance_transform_edt(out!=94);d95=distance_transform_edt(out!=95);out[hole&(d94<=d95)]=94;out[hole&(d95<d94)]=95
dpa=distance_transform_edt(~(L==66))*1.5
out[(out==95)&(dpa<14)]=98
out[Hm&(out==0)&(V>-20)]=97
mm=lambda m:distance_transform_edt(~m)*1.5
dLV,dRV,dRA,dLA,dRVOT,dPA=mm(out==4),mm(out==95),mm(out==94),mm(out==96),mm(out==98),mm(L==66)
out[(out==97)&(dLV<=9)&(dRV<=9)]=99                                   # interventricular septum
out[((out==94)&(dRV<=2.5))|((out==95)&(dRA<=2.5))]=100                 # tricuspid valve plane (approx.)
out[((out==98)&(dPA<=2.5))]=101                                        # pulmonary valve plane (approx.)
out[((out==96)&(dLV<=2.5))|((out==4)&(dLA<=2.5))]=102                  # mitral valve plane (approx.)
for i,n in ((4,'LV'),(96,'LA'),(94,'RA'),(95,'RV'),(98,'RVOT'),(97,'myo'),(99,'IVS'),(100,'TV'),(101,'PV'),(102,'MV')):
    m=out==i;print(n,int(m.sum()),'HU median',int(np.median(v[sl][m])) if m.any() else '-')
new=lab.copy();sub=new[sl];sub[Hm]=0;sub[out>0]=out[out>0];new[sl]=sub
np.save('labels_chambers.npy',new)
