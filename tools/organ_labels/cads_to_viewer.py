# Map the CADS output back into the viewer volume (z head->feet, y anterior->posterior, x right->left, same crop) and compare with the previous chamber labels.
import sys,numpy as np,nibabel as nib
from PIL import Image
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
lic=sys.argv[1] if len(sys.argv)>1 else 'open'
k0=262
c=np.asanyarray(nib.load(f'{S}/cads_out/chest_{lic}.nii.gz').dataobj).astype(np.uint8)
full=np.zeros((255,255,517),np.uint8);full[:,:,k0:k0+c.shape[2]]=c
a=np.asanyarray(nib.load('/home/claude/seg/ct1349.nii.gz').dataobj)
def tov(x): return x.transpose(2,1,0)[::-1,::-1,::-1]
im=tov(a);body=im>-400;ys=np.where(body.any(axis=(0,2)))[0];xs=np.where(body.any(axis=(0,1)))[0]
y0,y1=max(ys.min()-4,0),min(ys.max()+5,255);x0,x1=max(xs.min()-4,0),min(xs.max()+5,255)
cv=np.ascontiguousarray(tov(full)[:,y0:y1,x0:x1]);V=np.load(f'{S}/v.npy');assert cv.shape==V.shape,(cv.shape,V.shape)
np.save(f'{S}/cads_{lic}_view.npy',cv)
L=np.load('/home/claude/seg/labels_chambers.npy')
cur={'LV':[4,102],'LA':[96],'RA':[94,100],'RV':[95,98,101],'myo':[97,99],'PA':[66]}
cad={'LV':[5],'LA':[4],'RA':[6],'RV':[7],'myo':[3],'PA':[8]}
for k in cur:
    A=np.isin(L,cur[k]);B=np.isin(cv,cad[k]);print(f'{k:4s} current {A.sum()*3.375/1000:7.1f} mL  CADS {B.sum()*3.375/1000:7.1f} mL  Dice {2*(A&B).sum()/max(A.sum()+B.sum(),1):.2f}')
# where do current RA/RV voxels go in CADS?
for k in ['RA','RV']:
    A=np.isin(L,cur[k]);u,n=np.unique(cv[A],return_counts=True);print(k,'->',dict(zip(u.tolist(),(n/A.sum()).round(2).tolist())))
