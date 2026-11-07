import numpy as np,sys
from scipy.ndimage import map_coordinates
from PIL import Image,ImageDraw
S='.';name=sys.argv[1];out=sys.argv[2]
v=np.load(name).astype(np.float32);m=np.load('hyb_meta.npy');O,c,A=m[:3],m[3:6],m[6:]
c_loc=c-O
ant=np.array([0,-1.,0]);u1=ant-ant.dot(A)*A;u1/=np.linalg.norm(u1);u2=np.cross(A,u1)
lft=np.array([1.,0,0]);u2=lft-lft.dot(A)*A-(lft-lft.dot(A)*A).dot(u1)*u1;u2/=np.linalg.norm(u2)
def plane(cc,u,w,fov=70,Nn=280):
    s=np.linspace(-fov/2,fov/2,Nn);Ss,Tt=np.meshgrid(s,s);P=cc+Ss[...,None]*u+Tt[...,None]*w
    return map_coordinates(v,[P[...,2]/.5-.5,P[...,1]/.5-.5,P[...,0]/.5-.5],order=1,cval=-1000)
win=lambda a,l=250,w=900:np.clip((a-(l-w/2))/w*255,0,255).astype(np.uint8)
def tile(img,txt):
    im=Image.fromarray(win(img)).convert('RGB');ImageDraw.Draw(im).text((4,4),txt,fill=(255,200,60));return np.array(im)
t=[]
for dz in (-12,-4,0,6,14):
    t.append(tile(plane(c_loc+np.array([0,0,-1.])*dz,np.array([1.,0,0]),np.array([0,1.,0]),fov=100),f'axial {dz:+}'))
for hh in (-3,0,5,10,15):
    t.append(tile(plane(c_loc+A*hh,u2,-u1,fov=60),f'SA h={hh}'))
for ang in (0,40,80,120,160):
    r=np.radians(ang);d=np.cos(r)*u1+np.sin(r)*u2
    t.append(tile(plane(c_loc+A*10,d,-A,fov=80),f'LA {ang}'))
rows=[np.concatenate(t[i:i+5],1) for i in range(0,len(t),5)]
Image.fromarray(np.concatenate(rows,0)).save(out)
np.save('hyb_frame.npy',np.stack([u1,u2,A,c_loc]))
