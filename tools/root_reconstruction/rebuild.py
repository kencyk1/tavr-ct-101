import numpy as np,json,os
from scipy.ndimage import map_coordinates,gaussian_filter,gaussian_filter1d
sh=np.load('hyb_sh.npy').astype(np.float32);v3=np.load('hyb_vol3.npy').astype(np.float32);lab=np.load('hyb_lab4.npy')
u1,u2,A,c=np.load('hyb_frame.npy');R=np.load('hyb_Rfin.npy');HS=np.arange(-6,40.25,0.25);T=json.load(open('hyb_truth.json'))
rt=np.load('routes.npz');N=220;VS=.5;lum=T['lumenHU']
up=np.load('hyb_up.npy').astype(np.float32)
HS2=np.arange(40.25,70.25,0.25);TH_=np.radians(np.arange(0,360,1.0));RS=np.arange(0,30,0.2)
dd=np.cos(TH_)[:,None,None]*u1+np.sin(TH_)[:,None,None]*u2;Rx=np.zeros((len(HS2),360))
for i,h in enumerate(HS2):
    P=c+A*h+dd*RS[None,:,None];val=map_coordinates(up,[P[...,2]/.5-.5,P[...,1]/.5-.5,P[...,0]/.5-.5],order=1,cval=-1000)
    b=val<250;b[:,0]=False;r=np.where(b.any(1),RS[b.argmax(1)],RS[-1]);Rx[i]=np.minimum(r,np.percentile(r,50)*1.12)
from scipy.ndimage import median_filter
Rx=median_filter(Rx,(5,15),mode='wrap');Rx=gaussian_filter1d(gaussian_filter1d(Rx,4,axis=1,mode='wrap'),3,axis=0)
R=np.vstack([R,Rx]);HS=np.concatenate([HS,HS2]);del up
# ---- root surface Rs(theta,h): real lumen, extra smoothing; basal ring = fitted ellipse ----
Rs=gaussian_filter1d(gaussian_filter1d(R,3,axis=1,mode='wrap'),4,axis=0,mode='nearest')
cc=np.array(T['c']);off=cc-c;ox,oy=off@u1,off@u2
sa,sb=T['max']/2,T['min']/2
# ellipse major-axis direction: recover from hyb8 fit (store); recompute quickly from Rs near h=0
TH=np.radians(np.arange(0,360,1.0))
i0=int(np.argmin(abs(HS-0)));r0=R[[int(np.argmin(abs(HS-h))) for h in np.arange(-2,1.01,0.25)]].mean(0)
F=np.fft.rfft(r0);F[7:]=0;r0=np.fft.irfft(F,360);x=r0*np.cos(TH)-ox;y=r0*np.sin(TH)-oy
C=np.cov(np.stack([x,y]));w,V=np.linalg.eigh(C);ea=V[:,1];phi=np.arctan2(ea[1],ea[0])
# ray from c at angle th hits ellipse centred (ox,oy): solve |M(t*d - o)|=1
def ell_r(th):
    d=np.stack([np.cos(th),np.sin(th)]);cph,sph=np.cos(phi),np.sin(phi)
    def tr(v):return np.stack([(v[0]*cph+v[1]*sph)/sa,(-v[0]*sph+v[1]*cph)/sb])
    D=tr(d);O=tr(np.array([ox,oy])[:,None]*np.ones_like(th))
    a=(D**2).sum(0);b=-2*(D*O).sum(0);cq=(O**2).sum(0)-1
    return (-b+np.sqrt(b*b-4*a*cq))/(2*a)
re_=ell_r(TH)
for i,h in enumerate(HS):
    wgt=np.clip((3.0-h)/2.5,0,1)        # ellipse below 0.5 mm, blend to real by 3 mm
    sc=1-0.03*np.clip(-h/6,0,1)
    Rs[i]=wgt*re_*sc+(1-wgt)*Rs[i]
np.save('Rs.npy',Rs)
ax=((np.arange(N)+.5)*VS).astype(np.float32);Z,Y,X=np.meshgrid(ax,ax,ax,indexing='ij')
dX,dY,dZ=X-c[0],Y-c[1],Z-c[2];H=dX*A[0]+dY*A[1]+dZ*A[2];Q1=dX*u1[0]+dY*u1[1]+dZ*u1[2];Q2=dX*u2[0]+dY*u2[1]+dZ*u2[2]
Rr=np.hypot(Q1,Q2);THd=np.degrees(np.arctan2(Q2,Q1))%360;del dX,dY,dZ,Q1,Q2
RsV=map_coordinates(Rs,[np.clip((H-HS[0])/0.25,0,len(HS)-1),THd],order=1,mode='grid-wrap').astype(np.float32)
d_root=(Rr-RsV).astype(np.float32);del THd,RsV
wh=np.clip(np.minimum((H+6)/2.0,(70-H)/6.0),0,1).astype(np.float32)
d_root=np.where(wh>0,d_root,99).astype(np.float32)
# ---- coronary tubes ----
def tubes(paths,X,Y,Z):
    dc=np.full(X.shape,99,np.float32);lc=np.zeros(X.shape,np.uint8)
    for P,r0_,r1_,labv in paths:
        P=np.asarray(P,float);seg=np.linalg.norm(np.diff(P,axis=0),axis=1);s=np.concatenate([[0],np.cumsum(seg)]);L=s[-1]
        rad=lambda t:r0_+(r1_-r0_)*np.clip(t/L,0,1)
        for k in range(len(P)-1):
            a,b=P[k],P[k+1];rmax=max(rad(s[k]),rad(s[k+1]))+3
            lo=np.maximum((np.minimum(a,b)-rmax)/VS,0).astype(int);hi=np.minimum((np.maximum(a,b)+rmax)/VS+1,N).astype(int)
            sl=(slice(lo[2],hi[2]),slice(lo[1],hi[1]),slice(lo[0],hi[0]))
            v=b-a;wx,wy,wz=X[sl]-a[0],Y[sl]-a[1],Z[sl]-a[2];t=np.clip((wx*v[0]+wy*v[1]+wz*v[2])/max(v@v,1e-9),0,1)
            dist=np.sqrt((wx-t*v[0])**2+(wy-t*v[1])**2+(wz-t*v[2])**2)-rad(s[k]+t*seg[k])
            better=dist<dc[sl];dc[sl]=np.where(better,dist,dc[sl]);lc[sl]=np.where(better,labv,lc[sl])
    return dc,lc
lm_in=rt['lm_in'];rca_in=rt['rca_in']
LAD=np.vstack([lm_in,rt['LAD']]);LCX=np.vstack([lm_in,rt['LCX']]);RCA=np.vstack([rca_in,rt['RCA']])
d_cor,l_cor=tubes([(LAD,2.1,1.5,8),(LCX,2.1,1.4,8),(RCA,1.9,1.5,9)],X,Y,Z)
# ---- composite ----
lp=gaussian_filter(sh,2.0);tex=sh-lp
du=np.minimum(d_root,d_cor);near_c=d_cor<d_root
fw=np.clip((H+3.0)/3.0,0,1).astype(np.float32)   # wall fades out below the annulus (LVOT)
tw=np.where(near_c,0.7,1.2).astype(np.float32)
s=lambda x,w:np.clip(0.5-x/w,0,1)
a_in=s(du,0.5);a_w=np.clip(s(du-tw,0.5)-a_in,0,1);a_o=np.clip(s(du-tw-1.4,0.9)-a_in-a_w,0,1)
fwc=np.where(near_c,1.0,fw);a_w*=fwc;a_o*=fwc
fl=np.clip((H+6.0)/4.0,0,1)*wh
noise=sh-gaussian_filter(sh,1.2)
lv_real=gaussian_filter(sh,2.5)
root_in=(lum+1.7*noise)*fl+(lv_real+noise)*(1-fl)
cor_in=lum+1.5*noise
Lin=np.where(d_root<0,root_in,cor_in)
Wv=np.where(near_c,55,125)+0.6*tex
Ov=sh-0.85*np.clip(lp-np.where(near_c,-40,140),0,None)
rest=1-a_in-a_w-a_o
base=sh
gate=np.where(near_c,1.0,wh)
a_in*=gate;a_w*=gate;a_o*=gate;rest=1-a_in-a_w-a_o
out=a_in*Lin+a_w*Wv+a_o*Ov+rest*base
# keep real & leaflet calcium
out=np.where((sh>650)&(d_root>-0.5),np.maximum(out,sh),out)   # keep real wall / coronary calcium
# ---- repaint cusps fitted to the new wall, and leaflet calcium ----
Hc,Hf,TL,HU_LEAF=20.0,16.0,2.2,95.0
cen={'RCC':-23.0,'LCC':118.0,'NCC':-105.0}
comm={('RCC','LCC'):47.5,('LCC','NCC'):186.5,('NCC','RCC'):296.0}
bounds={'RCC':(comm[('NCC','RCC')]-360,comm[('RCC','LCC')]),'LCC':(comm[('RCC','LCC')],comm[('LCC','NCC')]),'NCC':(comm[('LCC','NCC')]-360,comm[('NCC','RCC')]-360)}
def Rlook(h,thd):return map_coordinates(Rs,[np.clip((h-HS[0])/0.25,0,len(HS)-1),(thd%360)],order=1,mode='grid-wrap')
dX,dY,dZ=X-c[0],Y-c[1],Z-c[2];Q1=dX*u1[0]+dY*u1[1]+dZ*u1[2];Q2=dX*u2[0]+dY*u2[1]+dZ*u2[2];THd=np.degrees(np.arctan2(Q2,Q1))%360;del dX,dY,dZ,Q1,Q2
roi=(H>-1.5)&(H<Hc+1)&(Rr<27);idx=np.nonzero(roi);h_=H[idx];th_=THd[idx];r_=Rr[idx];sd=np.full(h_.shape,99.0,np.float32);leaf_par={}
for nm in ('RCC','LCC','NCC'):
    lo_,hi_=bounds[nm];ctr=cen[nm];dd=((th_-ctr+180)%360)-180;inside=(dd>=lo_-ctr)&(dd<=hi_-ctr)
    half=np.where(dd<0,ctr-lo_,hi_-ctr);tt=np.clip(np.abs(dd)/half,0,1);hatt=Hc*(1-np.cos(np.pi/2*tt));ok=inside&(hatt<Hf-1.0)
    ss=(h_-hatt)/np.maximum(Hf-hatt,1.0);ok&=(ss>=0)&(ss<=1);ratt=Rlook(hatt,th_);rfree=4.5+1.0*np.cos(np.pi*tt)
    sc=np.clip(ss,0,1)**0.75;rl=ratt*(1-sc)+rfree*sc;sd=np.where(ok,np.minimum(sd,np.abs(r_-rl)-TL/2),sd);leaf_par[nm]=(lo_,hi_,ctr)
aL=np.zeros(out.shape,np.float32);aL[idx]=np.clip(0.5-sd/0.4,0,1);aL=np.clip(gaussian_filter(aL,0.85)*1.15,0,1)
out=out.astype(np.float32)*(1-aL)+(HU_LEAF+0.5*noise)*aL
rng=np.random.default_rng(5);aC=np.zeros(out.shape,np.float32);hC=np.zeros(out.shape,np.float32)
def leaf_pt(nm,fa,sv):
    lo_,hi_,ctr=leaf_par[nm];half=(ctr-lo_) if fa<0 else (hi_-ctr);thd=ctr+fa*half;tt=abs(fa);hatt=Hc*(1-np.cos(np.pi/2*tt));hh=hatt+sv*max(Hf-hatt,1)
    ra=float(Rlook(np.array([hatt]),np.array([thd]))[0]);rf=4.5+np.cos(np.pi*tt);rr_=ra*(1-sv**.75)+rf*sv**.75;t_=np.radians(thd);return c+A*hh+rr_*(np.cos(t_)*u1+np.sin(t_)*u2)
nod=[]
for nm,nn in (('NCC',3),('LCC',2),('RCC',2)):
    for _ in range(nn):
        p=leaf_pt(nm,rng.uniform(-.6,.6),rng.uniform(.15,.7));r0_=rng.uniform(1.0,1.8);nod.append((p,r0_))
        for _ in range(2):nod.append((p+rng.standard_normal(3)/np.sqrt(3)*r0_*.8,r0_*rng.uniform(.4,.7)))
for p,r0_ in nod:
    lo_=np.maximum(((p-r0_-2)/VS).astype(int),0);hi_=np.minimum(((p+r0_+2)/VS).astype(int)+1,N);sl=(slice(lo_[2],hi_[2]),slice(lo_[1],hi_[1]),slice(lo_[0],hi_[0]))
    k=np.sqrt((X[sl]-p[0])**2+(Y[sl]-p[1])**2+(Z[sl]-p[2])**2)/r0_;al=np.clip(0.5-(k-1)*r0_/0.35,0,1);huv=rng.uniform(900,1200)
    aC[sl]=np.maximum(aC[sl],al);hC[sl]=np.where(al>0,huv,hC[sl])
aC=np.clip(gaussian_filter(aC,0.9)*1.2,0,1);hC=np.maximum(gaussian_filter(hC,0.9)/np.maximum(gaussian_filter((hC>0).astype(np.float32),0.9),1e-3),0)
out=np.maximum(out,out*(1-aC)+hC*aC)
# labels
lab2=lab.copy();lab2[np.isin(lab2,[3,8,9,10])]=0;lab2[aL>0.35]=3;lab2[aC>0.4]=10
lab2[(d_cor<0)&(d_root>0)]=l_cor[(d_cor<0)&(d_root>0)]
out=np.clip(out,-1024,3071).astype(np.int16);np.save('reb_vol.npy',out);np.save('reb_lab.npy',lab2);np.save('d_root.npy',d_root);np.save('d_cor.npy',d_cor)
# ---- truth: ostium lower edges (centerline crossing the root surface) ----
def ostium(P,r0_):
    rel=P-c;h=rel@A;th=np.degrees(np.arctan2(rel@u2,rel@u1))%360;rr=np.hypot(rel@u1,rel@u2)
    rs=map_coordinates(Rs,[np.clip((h-HS[0])/0.25,0,len(HS)-1),th],order=1,mode='grid-wrap');dd=rr-rs
    k=int(np.argmax(dd>0));f=-dd[k-1]/(dd[k]-dd[k-1]);hc=h[k-1]+f*(h[k]-h[k-1]);thc=th[k-1]+f*(((th[k]-th[k-1]+180)%360)-180)
    return round(float(hc-r0_),1),round(float(thc if thc<180 else thc-360),0)
lm=ostium(LAD,2.1);rc=ostium(RCA,1.9)
# STJ & sinus diameters from Rs
iH=lambda h:int(np.argmin(abs(HS-h)))
sinus=[2*Rs[iH(h)].mean() for h in np.arange(6,18,0.5)];hs=np.arange(6,18,0.5)[int(np.argmax(sinus))]
stj_h=None
for h in np.arange(hs,36,0.25):
    if Rs[iH(h)].mean()<=Rs[iH(h+0.25)].mean() and h>hs+6: stj_h=h;break
T['coronaries']={'LM':{'lower_edge_mm':lm[0],'angle_deg':lm[1]},'RCA':{'lower_edge_mm':rc[0],'angle_deg':rc[1]}}
T['sinus']=dict(h=float(hs),mean_d=round(float(max(sinus)),1));T['stj']=dict(h=float(stj_h),mean_d=round(float(2*Rs[iH(stj_h)].mean()),1),height_above_annulus=float(stj_h))
json.dump(T,open('reb_truth.json','w'));print(json.dumps({k:T[k] for k in ('coronaries','sinus','stj','max','min','area','per')}))
