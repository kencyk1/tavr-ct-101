import numpy as np
from scipy.ndimage import map_coordinates,gaussian_filter
from PIL import Image,ImageDraw
v=np.load('v.npy').astype(np.float32)
c_ann=np.array([207.0,142.3,255.6]);n_lv=np.array([0.838,0.134,0.528]);n_lv/=np.linalg.norm(n_lv);A=-n_lv
VS=.5;N=220;Cout=c_ann+A*6;O=Cout-N*VS/2
ax=((np.arange(N)+.5)*VS).astype(np.float32)
Z,Y,X=np.meshgrid(ax+O[2],ax+O[1],ax+O[0],indexing='ij')
up=map_coordinates(v,[Z/1.5-.5,Y/1.5-.5,X/1.5-.5],order=3,mode='nearest').astype(np.float32)
del Z,Y,X
sh=up+0.9*(up-gaussian_filter(up,3.0))
np.save('hyb_up.npy',up.astype(np.int16));np.save('hyb_sh.npy',np.clip(sh,-1024,3071).astype(np.int16))
np.save('hyb_meta.npy',np.concatenate([O,c_ann,A]))
print('O',O,'stats lumen',np.percentile(up[(up>250)&(up<700)],[10,50,90]))
