# Re-sample the (CADS-based) chamber labels into the existing root reconstruction without rebuilding its CT.
import numpy as np,json,os
from PIL import Image
from scipy.ndimage import map_coordinates
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad';OUT='/home/claude/ctp/data'
T=json.load(open(S+'/root_v3_truth.json'));O=np.array(T['cube_origin_world']);N=T['shape'][0];VS=T['voxel_mm']
labC=np.load('/home/claude/seg/labels_chambers.npy');old=np.load(S+'/root_v3_lab.npy');vol=np.load(S+'/root_v3_vol.npy')
ax=((np.arange(N)+.5)*VS).astype(np.float32);Z,Y,X=np.meshgrid(ax,ax,ax,indexing='ij')
labU=map_coordinates(labC,[(Z+O[2])/1.5-.5,(Y+O[1])/1.5-.5,(X+O[0])/1.5-.5],order=0,mode='nearest').astype(np.uint8)
labU[np.isin(labU,[1,2,3,5,7])]=0
keep=np.isin(old,[1,2,3,5,7,8,9,10])          # lumen, annulus, LVOT, cusps, calcium, coronaries built by root_v3
lab=np.where(keep,old,labU).astype(np.uint8)
print('changed voxels',int((lab!=old).sum()))
np.save(S+'/root_v3_lab.npy',lab)
u=(vol.astype(np.int32)+1024).astype(np.uint16);files=[]
for i,z0 in enumerate(range(0,N,64)):
    z1=min(z0+64,N);cu=u[z0:z1].reshape(-1,N);ll=lab[z0:z1].reshape(-1,N)
    rgb=np.stack([(cu>>8).astype(np.uint8),(cu&255).astype(np.uint8),ll],-1);fn=f'{OUT}/root_{i:02d}.png'
    Image.fromarray(rgb,'RGB').save(fn,optimize=True,compress_level=9);files.append([z0,z1,os.path.getsize(fn)])
print(json.dumps(files))
