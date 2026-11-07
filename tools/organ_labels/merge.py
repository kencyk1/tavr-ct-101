import nibabel as nib,numpy as np,json,sys,gzip,struct
from totalsegmentator.map_to_binary import class_map
from scipy.ndimage import binary_opening
cm=class_map['total'];segf=sys.argv[1]
S='/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
# beginner groups: (id, name, short, color, [TotalSegmentator class names], showName)
G=[]
def g(i,n,s,c,names,show=True):G.append((i,n,s,c,names,show))
art=[220,70,70];vein=[80,120,230]
g(20,'Right lung','R lung',[110,170,230],['lung_upper_lobe_right','lung_middle_lobe_right','lung_lower_lobe_right'])
g(21,'Left lung','L lung',[110,170,230],['lung_upper_lobe_left','lung_lower_lobe_left'])
g(22,'Heart','Heart',[230,120,150],['heart'])
g(23,'Trachea','Trachea',[150,200,220],['trachea'])
g(24,'Esophagus','Esophagus',[230,170,120],['esophagus'])
g(25,'Thyroid gland','Thyroid',[240,140,90],['thyroid_gland'])
g(26,'Liver','Liver',[180,110,70],['liver'])
g(27,'Gallbladder','Gallbladder',[120,190,90],['gallbladder'])
g(28,'Spleen','Spleen',[170,100,200],['spleen'])
g(29,'Stomach','Stomach',[230,190,110],['stomach'])
g(30,'Pancreas','Pancreas',[240,210,120],['pancreas'])
g(31,'Right kidney','R kidney',[230,150,80],['kidney_right'])
g(32,'Left kidney','L kidney',[230,150,80],['kidney_left'])
g(33,'Adrenal gland','Adrenal',[240,200,60],['adrenal_gland_right','adrenal_gland_left'],False)
g(34,'Small bowel','Small bowel',[220,200,140],['small_bowel'])
g(35,'Duodenum','Duodenum',[220,200,140],['duodenum'])
g(36,'Colon','Colon',[200,170,120],['colon'])
g(37,'Urinary bladder','Bladder',[240,230,120],['urinary_bladder'])
g(38,'Prostate','Prostate',[210,150,170],['prostate'])
g(39,'Spinal cord','Cord',[250,240,180],['spinal_cord'])
g(40,'Superior vena cava','SVC',vein,['superior_vena_cava'])
g(41,'Inferior vena cava','IVC',vein,['inferior_vena_cava'])
g(42,'Brachiocephalic (innominate) artery','Innominate',art,['brachiocephalic_trunk'])
g(43,'Right subclavian artery','R subclavian',art,['subclavian_artery_right'])
g(44,'Left subclavian artery','L subclavian',art,['subclavian_artery_left'])
g(45,'Right common carotid artery','R carotid',art,['common_carotid_artery_right'])
g(46,'Left common carotid artery','L carotid',art,['common_carotid_artery_left'])
g(47,'Brachiocephalic vein','BCV',vein,['brachiocephalic_vein_left','brachiocephalic_vein_right'],False)
g(48,'Pulmonary vein','Pulm vein',vein,['pulmonary_vein'],False)
g(49,'Left atrial appendage','LAA',[230,120,150],['atrial_appendage_left'])
g(50,'Portal / splenic vein','Portal v',vein,['portal_vein_and_splenic_vein'],False)
g(51,'Right common/external iliac vein','R iliac v',vein,['iliac_vena_right'],False)
g(52,'Left common/external iliac vein','L iliac v',vein,['iliac_vena_left'],False)
g(53,'Sternum','Sternum',[235,225,200],['sternum'])
g(54,'Rib','Rib',[235,225,200],[f'rib_{s}_{i}' for s in ('left','right') for i in range(1,13)],False)
g(55,'Costal cartilage','Cartilage',[220,215,190],['costal_cartilages'],False)
g(56,'Clavicle','Clavicle',[235,225,200],['clavicula_left','clavicula_right'])
g(57,'Scapula','Scapula',[235,225,200],['scapula_left','scapula_right'])
g(58,'Humerus','Humerus',[235,225,200],['humerus_left','humerus_right'])
g(59,'Hip bone (ilium)','Ilium',[235,225,200],['hip_left','hip_right'])
g(60,'Femur','Femur',[235,225,200],['femur_left','femur_right'])
g(61,'Sacrum','Sacrum',[235,225,200],['sacrum','vertebrae_S1'])
g(62,'Gluteal muscles','Gluteal',[200,120,110],['gluteus_maximus_left','gluteus_maximus_right','gluteus_medius_left','gluteus_medius_right','gluteus_minimus_left','gluteus_minimus_right'])
g(63,'Paraspinal (erector spinae) muscle','Paraspinal',[200,120,110],['autochthon_left','autochthon_right'])
g(64,'Iliopsoas muscle','Psoas',[200,120,110],['iliopsoas_left','iliopsoas_right'])
g(65,'Skull','Skull',[235,225,200],['skull'],False)
vid=70
for lv in ['C1','C2','C3','C4','C5','C6','C7']+[f'T{i}' for i in range(1,13)]+[f'L{i}' for i in range(1,6)]:
    g(vid,f'{lv} vertebra',lv,[235,225,200],[f'vertebrae_{lv}']);vid+=1
inv={v:k for k,v in cm.items()}
seg=np.asanyarray(nib.load(segf).dataobj).astype(np.uint8)          # (x,y,z) RAS
lut=np.zeros(256,np.uint8)
for i,n,s,c,names,show in G:
    for nm in names:
        if nm in inv: lut[inv[nm]]=i
# TotalSegmentator aorta / iliac arteries fill gaps in the dataset's TAVR labels
lut[inv['aorta']]=1;lut[inv['iliac_artery_left']]=6;lut[inv['iliac_artery_right']]=7
org=lut[seg].transpose(2,1,0)[::-1,::-1,::-1]
# same crop as the viewer volume
raw=gzip.open('/mnt/user-data/uploads/Downloads/tavr-ct-primer/source/ct1349.nii.gz').read();vo=int(struct.unpack('<f',raw[108:112])[0])
im=np.frombuffer(raw[vo:],np.int16,count=255*255*517).reshape(517,255,255)[::-1,::-1,::-1]
body=im>-400;ys=np.where(body.any(axis=(0,2)))[0];xs=np.where(body.any(axis=(0,1)))[0]
y0,y1=max(ys.min()-4,0),min(ys.max()+5,255);x0,x1=max(xs.min()-4,0),min(xs.max()+5,255)
org=org[:,y0:y1,x0:x1];print('organ map',org.shape,(y0,y1,x0,x1))
tav=np.load(S+'/lb.npy');assert tav.shape==org.shape
out=np.where(tav>0,tav,org).astype(np.uint8)
np.save('labels_merged.npy',out)
json.dump([dict(id=i,n=n,s=s,c=c,show=show) for i,n,s,c,names,show in G],open('groups.json','w'))
u,cnt=np.unique(out,return_counts=True);print(len(u),'labels in volume')
