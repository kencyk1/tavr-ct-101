# Run the CADS heart model (task 553) on the chest of case ct1349, on CPU.
# Weights: https://github.com/murong-xu/CADS/releases/download/cads-model-open_v1.0.0/Dataset553_Totalseg253.zip (CC BY-SA 4.0),
# unzipped under weights/open/. Usage: PS=128 python3 cads_infer.py open
# The trained 192^3 patch needs more than 6 GB RAM on CPU, so the patch is reduced to 128^3 (about 4 GB, 3.5 min on 2 cores).
# Paths point at the original build machine; edit S and the CT path to rerun.
import sys,os,time,numpy as np,nibabel as nib,torch
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
lic=sys.argv[1] if len(sys.argv)>1 else 'open'
im=nib.load('/home/claude/seg/ct1349.nii.gz');a=np.asanyarray(im.dataobj)
k0,k1=262,454                                   # chest: nifti k (feet->head); covers heart, PA, great vessels
crop=a[:,:,k0:k1];aff=im.affine.copy();aff[:3,3]=im.affine[:3,:3]@np.array([0,0,k0])+im.affine[:3,3]
fi=f'{S}/cads_in/chest_0000.nii.gz'
if not os.path.exists(fi): nib.save(nib.Nifti1Image(crop.astype(np.int16),aff),fi)
torch.set_num_threads(2)
os.environ['nnUNet_compile']='f'
from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
p=nnUNetPredictor(tile_step_size=0.5,use_gaussian=True,use_mirroring=False,perform_everything_on_device=False,device=torch.device('cpu'),verbose=False,verbose_preprocessing=False,allow_tqdm=True)
md=f'/home/claude/cads_eval/weights/{lic}/Dataset553_Totalseg253/nnUNetTrainerNoMirroring__nnUNetResEncUNetLPlans__3d_fullres'
p.initialize_from_trained_model_folder(md,'all','checkpoint_final.pth')
ps=int(os.environ.get('PS','128'));p.configuration_manager.configuration['patch_size']=[ps,ps,ps];print('patch',p.configuration_manager.patch_size,flush=True)
t=time.time()
from nnunetv2.imageio.simpleitk_reader_writer import SimpleITKIO
io=SimpleITKIO();img,props=io.read_images([fi])
seg=p.predict_single_npy_array(img,props,None,None,False)
io.write_seg(seg,f'{S}/cads_out/chest_{lic}.nii.gz',props)
print('done',time.time()-t,flush=True)
