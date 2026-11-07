import numpy as np,os,json
from PIL import Image
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
v=np.load(S+'/v.npy').astype(np.int32);lab=np.load('/home/claude/seg/labels_chambers.npy')
v[v<-900]=-1024;u=(v+1024).astype(np.uint16);Z,Y,X=u.shape;out=[]
for i,z0 in enumerate(range(0,Z,64)):
    z1=min(z0+64,Z);c=u[z0:z1].reshape(-1,X);l=lab[z0:z1].reshape(-1,X)
    rgb=np.stack([(c>>8).astype(np.uint8),(c&255).astype(np.uint8),l],-1);fn=f'/home/claude/ctp/data/ct_{i:02d}.png'
    Image.fromarray(rgb,'RGB').save(fn,optimize=True,compress_level=9);out.append([z0,z1,os.path.getsize(fn)])
print(json.dumps(out))
