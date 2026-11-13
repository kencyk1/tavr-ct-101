# usage: python3 cads_any.py <ct.nii.gz> <k0> <k1> <out.nii.gz>   (k = RAS z index range of the crop)
import sys,os,time,numpy as np,nibabel as nib,torch
src,k0,k1,out=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
im=nib.load(src);a=np.asanyarray(im.dataobj);crop=a[:,:,k0:k1];aff=im.affine.copy();aff[:3,3]=im.affine[:3,:3]@np.array([0,0,k0])+im.affine[:3,3]
fi=f'{S}/cads_in/crop_{os.path.basename(out)}_0000.nii.gz';os.makedirs(f'{S}/cads_in',exist_ok=True)
nib.save(nib.Nifti1Image(crop.astype(np.int16),aff),fi)
torch.set_num_threads(2);os.environ['nnUNet_compile']='f'
from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
from nnunetv2.imageio.simpleitk_reader_writer import SimpleITKIO
p=nnUNetPredictor(tile_step_size=0.5,use_gaussian=True,use_mirroring=False,perform_everything_on_device=False,device=torch.device('cpu'),verbose=False,verbose_preprocessing=False,allow_tqdm=True)
p.initialize_from_trained_model_folder('/home/claude/cads_eval/weights/open/Dataset553_Totalseg253/nnUNetTrainerNoMirroring__nnUNetResEncUNetLPlans__3d_fullres','all','checkpoint_final.pth')
p.configuration_manager.configuration['patch_size']=[128,128,128]
io=SimpleITKIO();img,props=io.read_images([fi]);t=time.time()
seg=p.predict_single_npy_array(img,props,None,None,False);io.write_seg(seg,out,props);print('done',time.time()-t,flush=True)
