# Heart chamber labels from the CADS open model (task 553: myocardium, LA, LV, RA, RV).
# Replaces the seeded-watershed + geometric RA/RV split of chambers2.py.
# Input: labels_merged.npy (TAVR dataset labels + TotalSegmentator organs), CADS segmentation in viewer orientation.
import sys,numpy as np
from scipy.ndimage import distance_transform_edt,binary_opening,label as cc
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
cadsf=sys.argv[1] if len(sys.argv)>1 else S+'/cads_open_view.npy'
v=np.load(S+'/v.npy').astype(np.float32);lab=np.load('labels_merged.npy');C=np.load(cadsf)
# CADS 553 ids -> page ids
MAP={5:4,4:96,6:94,7:95,3:97}            # LV, LA, RA, RV, LV myocardium
heartish=np.isin(lab,[0,4,22])          # CADS may only overwrite background, dataset LV and TotalSegmentator "heart"
tavr=np.isin(lab,[1,2,3,5])             # dataset aorta / root / valve / annulus always win
new=lab.copy()
new[lab==22]=0;new[(lab==4)&~np.isin(C,[5])&~tavr]=0
for k,i in MAP.items():
    m=(C==k)&heartish&~tavr
    # keep the largest connected piece only (drops stray islands)
    l,n=cc(m)
    if n>1:
        sz=np.bincount(l.ravel());sz[0]=0;m=l==sz.argmax()
    new[m]=i
# leftover TotalSegmentator heart tissue (RV free wall, atrial walls): muscle-density -> myocardium, fat -> background
left=(lab==22)&(new==0)
new[left&(v>-20)]=97
zz,yy,xx=np.where(np.isin(new,[4,94,95,96,97]));z0,z1,y0,y1,x0,x1=zz.min()-6,zz.max()+7,yy.min()-6,yy.max()+7,xx.min()-6,xx.max()+7
sl=(slice(z0,z1),slice(y0,y1),slice(x0,x1));o=new[sl]
mm=lambda m:distance_transform_edt(~m)*1.5
dpa=mm(lab[sl]==66)
o[(o==95)&(dpa<14)]=98                                                     # RVOT = RV within 14 mm of the pulmonary artery
dLV,dRV,dRA,dLA,dRVOT=mm(o==4),mm(o==95),mm(o==94),mm(o==96),mm(o==98)
o[(o==97)&(dLV<=9)&((dRV<=9)|(dRVOT<=9))]=99                             # interventricular septum
o[((o==94)&(dRV<=2.5))|((o==95)&(dRA<=2.5))]=100                          # tricuspid valve plane (RA-RV contact)
o[(o==98)&(dpa<=2.5)]=101                                                 # pulmonary valve plane (RVOT-PA contact)
o[((o==96)&(dLV<=2.5))|((o==4)&(dLA<=2.5))]=102                           # mitral valve plane (LA-LV contact)
new[sl]=o
for i,n in ((4,'LV'),(96,'LA'),(94,'RA'),(95,'RV'),(98,'RVOT'),(97,'myo'),(99,'IVS'),(100,'TV'),(101,'PV'),(102,'MV'),(22,'heart')):
    m=new==i;print(f'{n:5s} {m.sum()*3.375/1000:7.1f} mL',f'HU median {int(np.median(v[m]))}' if m.any() else '')
np.save('labels_chambers.npy',new)
