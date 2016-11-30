# TODO break into several modules in package


import os
import sys
import glob
import numpy as np
import nibabel as nbl
import scipy.io
import scipy.cluster
import warnings
import nibabel.freesurfer as fs
import nibabel.freesurfer.io as fsio
import matplotlib
matplotlib.use('Qt5Agg')
import pylab as pl

try:
    import gdist
except ImportError:
    warnings.warn('Geodesic distance module unavailable; please pip install gdist.')


SUBJECTS_DIR, SUBJECT, FREESURFER_HOME = [os.environ[key] for key in 'SUBJECTS_DIR SUBJECT FREESURFER_HOME'.split()]



#------------------------------Electrode files---------------------------------

def parse_asa_electrode_file(fname):
    "Parse an ASA electrode format file."
    contents = {'positions': [], 'labels': []}
    with open(fname, 'r') as fd:
        lines = (l for l in fd.readlines())
        # parse header
        for line in lines:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if line.startswith('ReferenceLabel'):
                contents['reference_label'] = parts[1]
            elif line.startswith('UnitPosition'):
                contents['unit_position'] = parts[1]
            elif line.startswith('NumberPositions'):
                contents['number_positions'] = int(parts[1])
            elif line.startswith('Positions'):
                break
            else:
                raise Exception('unknown header line: %r' % (line,))
        # parse positions
        for line, _ in zip(lines, range(contents['number_positions'])):
            contents['positions'].append(
                    [float(coord) for coord in line.strip().split()])
        # parse labels
        #assert next(lines).strip() == 'Labels'
        [contents['labels'].append(line.strip()) for line in lines]
    return contents


def test_parse_asa_electrode_file():
    contents = parse_asa_electrode_file('standard_1005.elc')
    assert contents['reference_label'] == 'avg'
    assert contents['unit_position'] == 'mm'
    assert contents['number_positions'] == 346
    assert len(contents['positions']) == 346
    assert contents['positions'][0] == [-86.0761, -19.9897, -47.9860]
    assert contents['positions'][-1] == [85.7939, -25.0093, -68.0310]
    assert len(contents['labels']) == 346
    assert contents['labels'][0] == 'LPA'
    assert contents['labels'][-1] == 'A2'


#-------------------------------Brain visa-------------------------------------
def write_brain_visa_surf(fname, v, f):
    vn = vertex_normals(v, f)
    with open(fname, 'w') as fd:
        fd.write('- %d\n' % len(vn))
        for (vx, vy, vz), (nx, ny, nz) in zip(v, vn):
            fd.write('%f %f %f %f %f %f\n' % (vx, vy, vz, nx, ny, nz))
        fd.write('- %d %d %d\n' % ((len(f),)*3))
        for i, j, k in f:
            fd.write('%d %d %d\n' % (i, j, k))

def convert_fs_to_brain_visa(fs_surf):
    v, f = fs.read_geometry(fs_surf)
    write_brain_visa_surf(fs_surf + '.tri', v, f)
    
 
#-------------------------------------BEM--------------------------------------
 
def convert_bem_to_tri():
    surfs_glob = '%s/%s/bem/watershed/*_surface-low' % (SUBJECTS_DIR, SUBJECT)
    for surf_name in glob.glob(surfs_glob):
        convert_fs_to_brain_visa(surf_name)


def gen_head_model():
    surfs_glob = '%s/%s/bem/watershed/*_surface-low.tri' % (SUBJECTS_DIR, SUBJECT)
    surfs = glob.glob(surfs_glob)

    if len(surfs) == 0:
        raise Exception('tri surfaces not found!')

    hm_base = '%s/%s/bem/head_model' % (SUBJECTS_DIR, SUBJECT)
    hm_temp = """# Domain Description 1.1

Interfaces 3

Interface Skull: "{0}_outer_skull_surface-low.tri"
Interface Cortex: "{0}_inner_skull_surface-low.tri"
Interface Head: "{0}_outer_skin_surface-low.tri"

Domains 4

Domain Scalp: Skull -Head
Domain Brain: -Cortex
Domain Air: Head
Domain Skull: Cortex -Skull
""".format('%s/%s/bem/watershed/%s' % (SUBJECTS_DIR, SUBJECT, SUBJECT))

    hm_geom = hm_base + '.geom'
    with open(hm_geom, 'w') as fd:
        fd.write(hm_temp)
    print ('%s written.' % (hm_geom,))

    hm_cond = hm_base + '.cond'
    with open(hm_cond, 'w') as fd:
        fd.write("""# Properties Description 1.0 (Conductivities)

Air         0.0
Scalp       1
Brain       1
Skull       0.03
""")
    print ('%s written.' % (hm_cond,))
    

            
#------------------------Freesurfer surface annotations------------------------

def annot_path(hemi, annot_name):
    annot_fname = '%s.%s.annot' % (hemi, annot_name)
    annot_path = os.path.join(SUBJECTS_DIR, SUBJECT, 'label', annot_fname)
    return annot_path

def read_annot(hemi, annot_name):
    return fs.read_annot(annot_path(hemi, annot_name))

def write_annot(hemi, annot_name, labels, ctab, names):
    return fs.write_annot(annot_path(hemi, annot_name), labels, ctab, names)

def read_surf(hemi, name):
    surf_fname = '%s.%s' % (hemi, name)
    surf_path = os.path.join(SUBJECTS_DIR, SUBJECT, 'surf', surf_fname)
    return fs.read_geometry(surf_path)

def read_lut(lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt'),key_mode='label'):
    from collections import OrderedDict
    f = open(lut_path,"r")
    l=list(f)
    f.close()
    ii=-1
    if key_mode=='label':
        names=OrderedDict()
        colors=OrderedDict()
        labels=[]
        for line in l:
            temp=line.split()
            try:
                label=int(temp[0])                
                ii+=1
                labels.append(label)
                names[labels[ii]]=temp[1]
                colors[labels[ii]]=[int(temp[2]),int(temp[3]),int(temp[4]),int(temp[5])]
            except:
                pass
    elif key_mode=='name':  
        labels=OrderedDict()
        colors=OrderedDict()
        names=[]
        for line in l:
            temp=line.split()
            try:
                label=int(temp[0])  
                ii+=1
                names.append(temp[1])
                labels[names[ii]]=label
                colors[names[ii]]=[int(temp[2]),int(temp[3]),int(temp[4]),int(temp[5])]
            except:
                pass
    return (labels,names,colors)

def annot_to_lut(hemi, annot_name, lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt')):
    _, ctab, names = read_annot(hemi, annot_name)
    with open(lut_path, 'w') as fd:
        for name, (r, g, b, a, id) in zip(names, ctab):
            fd.write('%d\t%s\t%d %d %d %d\n' % (id, name, r, g, b, a))

def rgb_to_fs_magic_number(rgb):
    return rgb[0]+256*rgb[1]+256*256*rgb[2]
    
def lut_to_annot_names_ctab(lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt'),labels=None):
    (_,names_dict,colors)=read_lut(lut_path=lut_path)
    if labels is None:
        labels=names_dict.keys()
    elif isinstance(labels, basestring):
        labels=np.array(labels.split()).astype('i').tolist()
    else:
        labels=np.array(labels).astype('i').tolist()
    names=[]
    ctab=[]
    for lbl in labels:
        names.append(names_dict[lbl])
        rgb=np.array(colors[lbl])[:3].astype('int64')
        magic_number=rgb_to_fs_magic_number(rgb)*np.ones((1,),dtype='int64')
        ctab.append(np.concatenate([rgb,np.zeros((1,),dtype='int64'),magic_number]))
    ctab=np.asarray(ctab).astype('int64')
    return (names,ctab)
    
def annot_names_to_labels(names,ctx=None,lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt')):
     (labels_dict,_,_)=read_lut(lut_path=lut_path,key_mode='name') 
     labels=[]
     if ctx=='lh' or ctx=='rh':
        ctx='ctx-'+ctx+'-'
     else:
        ctx=''
     for name in names:
         labels.append(labels_dict[ctx+name])
     return labels
    
def annot_to_conn_conf(hemi, annot_name, conn_conf_path):
    _, _, names = read_annot(hemi, annot_name)
    with open(conn_conf_path, 'w') as fd:
        for id, name in enumerate(names):
            fd.write('%d\t%s\n' % (id, name))

def read_input_labels(labels=None,hemi=None):
    if labels is not None:
        if isinstance(labels, basestring):
            #Set the target labels
            labels=np.array(labels.split()).astype('i').tolist() 
    else:
        labels=[]
    if hemi is not None:
        hemi=hemi.split()
        for h in hemi:
            if h=='lh':
                labels=labels+range(1000,1036)
            elif h=='rh':
                labels=labels+range(2000,2036)
    nLbl=len(labels)
    return (labels,nLbl)



#-----------------------------Freesurfer surfaces------------------------------

def tri_area(tri):
    i, j, k = np.transpose(tri, (1, 0, 2))
    ij = j - i
    ik = k - i
    return np.sqrt(np.sum(np.cross(ij, ik)**2, axis=1)) / 2.0  
    
def vertex_normals(v, f):
    vf = v[f]
    fn = np.cross(vf[:,1] - vf[:, 0], vf[:, 2] - vf[:, 0])
    vf = [set() for _ in range(len(v))]
    for i, fi in enumerate(f):
        for j in fi:
            vf[j].add(i)
    vn = np.zeros_like(v)
    for i, fi in enumerate(vf):
        fni = fn[list(fi)]
        norm = fni.sum(axis=0)
        norm /= np.sqrt((norm**2).sum())
        vn[i] = norm
    return vn

def compute_gdist_mat(surf_name='pial', max_distance=40.0):
    max_distance = float(max_distance) # in case passed from sys.argv
    for h in 'rl':
        surf_path = '%s/%s/surf/%sh.%s' % (SUBJECTS_DIR, SUBJECT, h, surf_name)
        v, f = fs.read_geometry(surf_path)
        mat_path = '%s/%s/surf/%sh.%s.gdist.mat' % (SUBJECTS_DIR, SUBJECT, h, surf_name)
        mat = gdist.local_gdist_matrix(v, f.astype('<i4'), max_distance=40.0)
        scipy.io.savemat(mat_path, {'gdist': mat})


def extract_subsurf(verts,faces,verts_mask):
    #These are the faces to keep...
    verts_out=verts[verts_mask,:]
    #These are the faces to keep...
    face_mask=np.c_[verts_mask[faces[:,0]],verts_mask[faces[:,1]],verts_mask[faces[:,2]]].all(axis=1)
    faces_out=faces[face_mask]
    #...but the old vertices' indexes of faces have to be transformed to the new verts_out_inds:
    verts_out_inds,=np.where(verts_mask)
    for iF in range(faces_out.shape[0]):
        for iV in range(3):
            faces_out[iF,iV],= np.where(faces_out[iF,iV]==verts_out_inds)
    return (verts_out, faces_out)  
    
def extract_mri_vol2subsurf(surf_path,annot_path,vol2surf_path,out_surf_path=None,out_annot_path=None,ctx=None,labels=None,lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt')):
   (verts, faces,volume_info) = fsio.read_geometry(surf_path,read_metadata=True)
   vol2surf=nbl.load(vol2surf_path)
   vol2surf=np.round(np.squeeze(vol2surf.get_data())).astype('i')
   lab, ctab, names = fs.read_annot(annot_path)
   if labels is None:
       labels=np.array(annot_names_to_labels(names,ctx=ctx,lut_path=lut_path))
   else:
       labels=np.array(labels.split()).astype('i')
   verts_mask=(v2s in labels for v2s in vol2surf)    
   (verts_out, faces_out)=extract_subsurf(verts,faces,verts_mask)
   if os.path.exists(str(out_surf_path)):
       fsio.read_geometry(out_surf_path,verts_out, faces_out,volume_info=volume_info)
   if os.path.exists(str(out_annot_path)):   
       lab=lab[verts_mask]
       fs.write_annot(out_annot_path,lab, ctab, names)
   return (verts_out, faces_out)     
              
#Concatenate surfaces of specific labels to create a single annotated surface
def aseg_surf_conc_annot(surf_path,out_surf_path,annot_path,labels,lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt')):
    names,ctab=lut_to_annot_names_ctab(lut_path=lut_path,labels=labels)
    labels=np.array(labels.split()).astype('i')
    out_verts=[]
    out_faces=[]
    lab=[]
    iL=-1
    nVerts=0
    names_out=[]
    ctab_out=[]
    for lbl in labels:
        l=int(lbl)
        this_surf_path=(surf_path+"-%06d" % (l))
        if os.path.exists(this_surf_path):
            indL,=np.where(labels==lbl)
            names_out.append(names[indL])
            ctab_out.append(ctab[indL,:])
            iL+=1
            (verts, faces,volume_info) = fsio.read_geometry(this_surf_path,read_metadata=True)
            faces=faces+nVerts #Update vertices indexes
            nVerts+=verts.shape[0]
            out_verts.append(verts)
            out_faces.append(faces)
            lab.append(iL*np.ones((verts.shape[0],),dtype='int64'))
    ctab_out=np.squeeze(np.array(ctab_out).astype('i'))        
    out_verts=np.vstack(out_verts)
    out_faces=np.vstack(out_faces)
    lab=np.hstack(lab)
    fsio.write_geometry(out_surf_path,out_verts,out_faces,create_stamp=None,volume_info=volume_info)  
    fs.write_annot(annot_path,lab, ctab_out, names_out)    


#It returns a sparse matrix of the connectivity among the vertices of a surface
#mode: "sparse" (default) or "2D"
def vertex_connectivity(v,f,mode="sparse",metric=None):
    #Get all pairs of vertex indexes that appear in each face
    f=np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]
    #Remove repetitions
    f=np.vstack(set(map(tuple,f)))
    #Mark all existing pairs to 1
    nV=v.shape[0]
    nF=f.shape[0]
    from scipy.sparse import csr_matrix
    if metric is None:
        con=csr_matrix((np.ones((nF,)),(f[:,0],f[:,1])),shape=(nV,nV))
        if mode!="sparse":
            #Create non-sparse matrix
            con=con.todense()
    else:
        from sklearn.metrics.pairwise import paired_distances
        d=paired_distances(v[f[:,0]],v[f[:,1]],metric)
        if mode=="sparse":
            #Create sparse matrix
            con=csr_matrix((d,(f[:,0],f[:,1])),shape=(nV,nV))  
    return con
    
    
#---------------------------------Volumes--------------------------------------  

#Separate the voxels of the outer surface of a structure, from the inner ones
#Default behavior: surface voxels retain their label, inner voxels get the label 0,
#and the input file is overwritten by the output
def vol_to_ext_surf_vol(in_vol_path, labels=None, hemi=None, out_vol_path=None,labels_surf=None,labels_inner='0'):
    #Set the target labels:
    (labels,nLbl)=read_input_labels(labels=labels,hemi=hemi)
    #Set the labels for the surfaces
    if labels_surf is None:
        labels_surf=labels
    else:
        #Read the surface labels
        labels_surf=np.array(labels_surf.split()).astype('i')
        #...and make sure there is one for each label
        if len(labels_surf)==1:
            labels_surf=np.repeat(labels_inner,nLbl).tolist()
        elif len(labels_surf)!=nLbl:
            print "Output labels for surface voxels are neither of length 1 nor of length equal to the one of target labels" 
        else:
            labels_surf=labels_surf.tolist()
    #Read the inner, non-surface labels
    labels_inner=np.array(labels_inner.split()).astype('i')
    #...and make sure there is one for each label
    if len(labels_inner)==1:
        labels_inner=np.repeat(labels_inner,nLbl).tolist()
    elif len(labels_inner)!=nLbl:
        print "Output labels for inner voxels are neither of length 1 nor of length equal to the one of the target labels"
    else:
        labels_inner=labels_inner.tolist()
    #Read the input volume...
    volume = nbl.load(in_vol_path)
    #...and get its data
    vol = volume.get_data()
    vol_shape=vol.shape
    #Neigbors' grid sharing a face
    borderGrid=np.c_[np.identity(3),-np.identity(3)].T.astype('i')
    nBorder=6 
    #Initialize output volume array
    out_vol=np.array(vol)
    #Initialize output indexes
    out_ijk=[]
    #For each target label:
    for iL in range(nLbl):
        #this label
        lbl=labels[iL]
        #Get the indexes of all voxels of this label:
        ii,jj,kk=np.where(vol==lbl)
        #and for each voxel
        for iV in range(ii.size):
            #indexes of this voxel:
            (i,j,k)=(ii[iV],jj[iV],kk[iV])
            #Create the neighbors' grid sharing a face
            ijk_grid=borderGrid+np.tile(np.array([i,j,k]),(nBorder,1))
            #Remove voxels outside the image
            inds_inside_image=np.all([(ijk_grid[:, 0]>=0), (ijk_grid[:, 0]<vol_shape[0]),
                                      (ijk_grid[:, 1]>=0), (ijk_grid[:, 1]<vol_shape[1]),
                                      (ijk_grid[:, 2]>=0), (ijk_grid[:, 2]<vol_shape[2])],
                                    axis=0)
            ijk_grid=ijk_grid[inds_inside_image,:]
            try:
                #If all face neighbors are of the same label...
                if np.all(vol[ijk_grid[:, 0],ijk_grid[:, 1],ijk_grid[:, 2]]==np.tile(vol[i,j,k],(nBorder,1))):
                    #...set this voxel to the corresponding inner target label
                    out_vol[i,j,k]=labels_inner[iL]
                else:
                    #...set this voxel to the corresponding surface target label
                    out_vol[i,j,k]=labels_surf[iL]
                    out_ijk.append([i,j,k])
            except ValueError: #empty grid
                print "Error at voxel ("+str(i)+","+str(j)+","+str(k)+") of label "+str(lbl)+":"
                print "It appears to have no common-face neighbors inside the image!"
                return     
    #Create the new volume and save it                                  
    out_volume=nbl.Nifti1Image(out_vol,volume.affine,header=volume.header)
    if out_vol_path==None:  
        #Overwrite volume
        out_vol_path=in_vol_path  
    #Save a new volume
    nbl.save(out_volume,out_vol_path) 
    #...and the output indexes that survived masking 
    out_ijk=np.vstack(out_ijk)
    filepath=os.path.splitext(out_vol_path)[0]
    np.save(filepath+"-idx.npy",out_ijk)
    np.savetxt(filepath+"-idx.txt",out_ijk,fmt='%d')
    
        
#Identify the voxels that our neighbors with a voxel distance vn, to a mask volume, 
#with a mask threshold of th
#Default behavior: we assume a binarized mask and set th=0.999,
#no neigbhors search, only looking at the exact voxel position, i.e., vn=0.        
#and mask voxels retain their label, no mask voxels get a label of 0        
def mask_to_vol(in_vol_path,mask_vol_path,out_vol_path=None,labels=None,hemi=None,vol2mask_path=None,vn=1,th=0.999,labels_mask=None,labels_nomask='0'):
    #Set the target labels:    
    (labels,nLbl)=read_input_labels(labels=labels,hemi=hemi)
    #Set the labels for the selected voxels
    if labels_mask is None:
        labels_mask=labels
    else:
        #Read the labels
        labels_mask=np.array(labels_mask.split()).astype('i')
        #...and make sure there is one for each label
        if len(labels_mask)==1:
            labels_mask=np.repeat(labels_mask,nLbl).tolist()
        elif len(labels_mask)!=nLbl:
            print "Output labels for selected voxels are neither of length 1 nor of length equal to the one of target labels" 
        else:
            labels_mask=labels_mask.tolist()
    #Read the excluded labels
    labels_nomask=np.array(labels_nomask.split()).astype('i')
    #...and make sure there is one for each label
    if len(labels_nomask)==1:
        labels_nomask=np.repeat(labels_nomask,nLbl).tolist()
    elif len(labels_nomask)!=nLbl:
        print "Output labels for excluded voxels are neither of length 1 nor of length equal to the one of the target labels"
    else:
        labels_nomask=labels_nomask.tolist()
    #Read the target volume...
    volume = nbl.load(in_vol_path)
    #...and get its data
    vol = volume.get_data()
     #...and its affine transform
    #ijk2xyz_vol = volume.affine
    #Read the mask volume...
    mask_vol = nbl.load(mask_vol_path)
    #...and get its data
    mask = mask_vol.get_data()
    mask_shape=mask.shape
    #...and invert its affine transform
    #xyz2ijk_mask = np.linalg.inv(mask_vol.affine)
    #Finally compute the transform from vol ijk to mask ijk:
    ijk2ijk=np.identity(4)
    #If vol and mask are not in the same space:
    if os.path.exists(str(vol2mask_path)):
        #read the xyz2xyz transform...
        xyz2xyz=np.loadtxt(vol2mask_path)
        #...and apply it to the inverse mask affine transform to get an ijk2ijk transform:
        ijk2ijk=volume.affine.dot(np.dot(xyz2xyz,np.linalg.inv(mask_vol.affine)))     
    #Construct a grid template of voxels +/- vn voxels around each ijk voxel,
    #sharing at least a corner
    grid=np.meshgrid(range(-vn,vn+1,1),range(-vn,vn+1,1),range(-vn,vn+1,1),indexing='ij')
    grid=np.c_[np.array(grid[0]).flatten(),np.array(grid[1]).flatten(),np.array(grid[2]).flatten()]
    nGrid=grid.shape[0]
    #Initialize the output volume
    out_vol=np.array(vol)
    #Initialize output indexes
    out_ijk=[]
    #For each target label:
    for iL in range(nLbl):
        lbl=labels[iL]   
        #Get the indexes of all voxels of this label:
        ii,jj,kk=np.where(vol==lbl)
        #and for each voxel
        for iV in range(ii.size):
            #indexes of this voxel:
            (i,j,k)=(ii[iV],jj[iV],kk[iV])
            #TODO if necessary: deal with voxels at the edge of the image, such as brain stem ones...
#           #Check if it is a border voxel:
#           if any([(i==0), (i==mask_shape[0]-1),(j==0), (j==mask_shape[0]-1),(k==0), (k==mask_shape[0]-1)]):
#               #set this voxel to the 0 label
#               mask_shape[i,j,k]=0 
#               continue           
            #...get the corresponding voxel in the mask volume:
            ijk=np.round(ijk2ijk.dot(np.array([i,j,k,1]))[:3]).astype('i')
            #Make sure this point is within image limits
            for cc in range(3):
                if ijk[cc]<0:
                    ijk[cc]=0
                elif ijk[cc]>=mask_shape[cc]: 
                    ijk[cc]=mask_shape[cc]-1
            #If this is a voxel to keep, set it so...
            if (mask[ijk[0],ijk[1],ijk[2]]>=th):
                out_vol[i,j,k]=labels_mask[iL]
                out_ijk.append([i,j,k])
            elif vn>0:
                #...if not, and as long as vn>0...
                #...check whether any of its vn neighbors is a mask voxel 
                #Generate the specific grid centered at the vertex ijk
                ijk_grid=grid+np.tile(ijk,(nGrid,1))
                #Remove voxels outside the mask volume
                indexes_within_limits=np.all([(ijk_grid[:, 0]>=0), (ijk_grid[:, 0]<mask_shape[0]),
                                              (ijk_grid[:, 1]>=0), (ijk_grid[:, 1]<mask_shape[1]),
                                              (ijk_grid[:, 2]>=0), (ijk_grid[:, 2]<mask_shape[2])],
                                            axis=0)
                ijk_grid=ijk_grid[indexes_within_limits,:]
                try:
                    #If none of these points is a mask point:
                    if (mask[ijk_grid[:, 0],ijk_grid[:, 1],ijk_grid[:, 2]]<th).all():
                        out_vol[i,j,k]=labels_nomask[iL]
                    else: #if any of them is a mask point:
                        out_vol[i,j,k]=labels_mask[iL]
                        out_ijk.append([i,j,k])
                except ValueError: #empty grid
                    print "Error at voxel ("+str(i)+","+str(j)+","+str(k)+"):"
                    print "It appears to have no common-face neighbors inside the image!"
                    return
            else:
                out_vol[i,j,k]=labels_nomask[iL]
    #Create the new volume and save it                                     
    out_volume=nbl.Nifti1Image(out_vol,volume.affine,header=volume.header)
    if out_vol_path==None:  
        #Overwrite volume
        out_vol_path=in_vol_path  
    nbl.save(out_volume,out_vol_path)
    #...and the output indexes that survived masking 
    out_ijk=np.vstack(out_ijk)
    filepath=os.path.splitext(out_vol_path)[0]
    np.save(filepath+"-idx.npy",out_ijk)
    np.savetxt(filepath+"-idx.txt",out_ijk,fmt='%d')
    
def label_with_dilation(to_label_nii_fname, dilated_nii_fname, out_nii_fname):
    "Label one nifti with its dilation, cf seeg-ct.sh"
    # TODO could make dilation with ndimage also.
    import nbl, scipy.ndimage
    mask = nbl.load(to_label_nii_fname)
    dil_mask = nbl.load(dilated_nii_fname)
    lab, n = scipy.ndimage.label(dil_mask.get_data())
    print('[label_with_dilation] %d objects found.' % (n, ))
    lab_mask_nii = nbl.nifti1.Nifti1Image(lab * mask.get_data(), mask.affine)
    nbl.save(lab_mask_nii, out_nii_fname)


def label_vol_from_tdi(tdi_nii_fname, out_fname, lo=0.5):
    "Make label volume from tckmap output."
    #Load tdi_ends volume:
    nii = nbl.load(tdi_nii_fname)
    #Copy its data...
    tdi = nii.get_data().copy()
    #and mask them to get the voxels of tract ends
    mask = tdi > lo
    #(all other voxels ->0)
    tdi[~mask] = 0
    #Assign them with integer labels starting from 1
    tdi[mask] = np.r_[1:mask.sum()+1]
    #Write tdi_lbl to file
    out_nii = nbl.nifti1.Nifti1Image(tdi, nii.affine)
    nbl.save(out_nii, out_fname)  
    
#It removes network nodes with zero connectivity, and returns a symmetric connectivity matrix
#Inputs:
#    - the tdi_lbl.nii volume path
#    -a .csv file path, output of Mrtrix3 tck2connectome
#Outputs:
#   - the symmetric no-zero-connections connectivity matrix saved as .npy
#   - the tdi_lbl.nii volume with the removed voxel nodes to 0 and the labels 
#    updated
#Optionally: if the tract length matrix is in the input, it is also processed 
def remove_zero_connectivity_nodes(node_vol_path,con_mat_path,tract_length_path=None):
    #Read input files:
    #Nodes' volume (nii):
    node_vol=nbl.load(node_vol_path)
    vol=node_vol.get_data()
    #Connectivity matrix (.csv)
    con=np.array(np.genfromtxt(con_mat_path,dtype='int64'))
    #Make it symmetric:
    con=con+con.T
    #Sum rows to get the total connectivity per node
    consum=np.sum(con,axis=0)
    #Index of nodes to keep
    ii=consum>0
    #Select only the specified columns and rows from the connectivity matrix:
    con=con[ii,:][:,ii]
    #Write output files
    np.save(os.path.splitext(con_mat_path)[0]+".npy", con)
    np.savetxt(con_mat_path, con)
    #If there is also a tract length file:
    if os.path.exists(str(tract_length_path)):
        #read it:
        con=np.array(np.genfromtxt(tract_length_path,dtype='int64'))
        #Select only the specified columns and rows from the connectivity matrix:
        con=con[ii,:][:,ii]
        #Write output files
        np.save(os.path.splitext(tract_length_path)[0]+".npy", con)
        np.savetxt(tract_length_path, con)
    else: 
        print tract_length_path+" is not a valid path"
    #Index of nodes to remove        
    ii,=np.where(~ii)
    ii=ii+1
    nKeep=con.shape[0]
    #Remove:
    for iR in ii:
        vol[vol==iR]=0
    #Update remaining indexes
    vol[vol>0]=np.r_[1:(nKeep+1)]
    #Write the updated volume file:
    node_vol=nbl.Nifti1Image(vol,node_vol.affine,header=node_vol.header)
    nbl.save(node_vol,node_vol_path)
    
#It receives a binary connectivity matrix, and outputs a node connectivity
#similarity or distance  matrix 
#con_mat_path: path to connectivity file
#metric: default "cosine"
#mode: "sim" or "dist" for similarity or distance output
def node_connectivity_metric(con_mat_path,metric="cosine", mode='sim', out_consim_path=None):
    from scipy.spatial.distance import pdist, squareform
    #Read iput file
    con=np.load(con_mat_path)
    #Calculate distance metric
    con=squareform(pdist(con, metric=metric))
    #If similarity is required,... 
    if mode=='sim':
        #...calculate it:
        con=1-con
    if out_consim_path is not None:
        np.save(out_consim_path,con)
    return con

    
def simple_label_config(aparc_fname, out_fname):
    "Rewrite label volume to have contiguous values like mrtrix' labelconfig."
    aparc = nbl.load(aparc_fname)
    vol = aparc.get_data()
    uval = np.unique(vol)
    uval_map = np.r_[:uval.max() + 1]
    uval_map[uval] = np.r_[:uval.size]
    uvol = uval_map[vol]
    uparc = nbl.nifti1.Nifti1Image(uvol, aparc.affine)
    nbl.save(uparc, out_fname)
  

 #-------------------------Surfaces from/to volumes----------------------------  

#Sample a volume of a specific label on a surface, by keeping  
#only those surface vertices, the nearest voxel of which is of the given label (+ of possibly additional target labels, such as white matter)
#Allow optionally for vertices within a given voxel distance vn from the target voxels
def sample_vol_on_surf(surf_path,vol_path,annot_path,out_surf_path,vox2rastkr_path, ctx=None,vn=1,add_lbl=[],lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt')):
    #Read the surfaces...
    (verts, faces,volume_info) = fsio.read_geometry(surf_path,read_metadata=True)
    #...and its annotation:
    lab, ctab, names = fs.read_annot(annot_path)
    #Get the region names of these labels:
    labels=annot_names_to_labels(names,ctx,lut_path)          
    nLbl=len(labels)           
    #Read the volume...
    volume = nbl.load(vol_path)
    #...and get its data
    vol = volume.get_data()
    vol_shape=vol.shape  
    #...and invert its vox2ras-tkr transform
    vox2rastkr=np.loadtxt(vox2rastkr_path)
    xyz2ijk = np.linalg.inv(vox2rastkr) 
    #Prepare grid if needed for possible use:
    if vn>0:
        grid=np.meshgrid(range(-vn,vn+1,1),range(-vn,vn+1,1),range(-vn,vn+1,1),indexing='ij')
        grid=np.c_[np.array(grid[0]).flatten(),np.array(grid[1]).flatten(),np.array(grid[2]).flatten()]
        nGrid=grid.shape[0]
    #Initialize the output mask:
    verts_out_mask=np.repeat(False,verts.shape[0])
    for iL in range(nLbl):
        if isinstance(ctx, basestring):
            print 'ctx-'+ctx+'-'+names[iL]
        else:
            print names[iL]
        #Form the target label by adding to the current label any additional labels, if any
        lbl=add_lbl+labels[iL]
        #Get the indexes of the vertices of this label:
        verts_lbl_inds,=np.where(lab[:]==iL)
        if verts_lbl_inds.size==0:
            continue
        #Apply the affine transform to the selected surface vertices, and get integer indices
        #of the corresponding nearest voxels in the vol
        ijk = np.round(xyz2ijk.dot(np.c_[verts[verts_lbl_inds,:], np.ones(verts_lbl_inds.size)].T)[:3].T).astype('i')
        #Get the labels of these voxels:
        surf_vxls=vol[ijk[:,0],ijk[:,1],ijk[:,2]]        
        #Vertex mask to keep: those that correspond to voxels of one of the target labels
        verts_keep,=np.where(np.in1d(surf_vxls,lbl)) #surf_vxls==lbl if only one target label
        verts_out_mask[verts_lbl_inds[verts_keep]]=True
        if vn>0:
            #These are now the remaining indexes to be checked for neighboring voxels
            verts_lbl_inds=np.delete(verts_lbl_inds,verts_keep)
            ijk=np.delete(ijk,verts_keep,axis=0)
            for iV in range(verts_lbl_inds.size):
                #Generate the specific grid centered at the voxel ijk
                ijk_grid=grid+np.tile(ijk[iV,:],(nGrid,1))
                #Remove voxels outside the volume
                indexes_within_limits=np.all([(ijk_grid[:, 0]>=0), (ijk_grid[:, 0]<vol_shape[0]),
                                              (ijk_grid[:, 1]>=0), (ijk_grid[:, 1]<vol_shape[1]),
                                              (ijk_grid[:, 2]>=0), (ijk_grid[:, 2]<vol_shape[2])],
                                             axis=0)
                ijk_grid=ijk_grid[indexes_within_limits,:]
                #Get the labels of these voxels:
                surf_vxls=vol[ijk_grid[:, 0],ijk_grid[:, 1],ijk_grid[:, 2]]
                #If any of the neighbors is of the target labels...
                if np.any(np.in1d(surf_vxls,lbl)): #surf_vxls==lbl if only one target label
                    #...include this vertex
                    verts_out_mask[verts_lbl_inds[iV]]=True
        #Vertex indexes to keep:
        verts_out_inds,=np.where(verts_out_mask)
        #These are the vertices to keep
        verts_out=verts[verts_out_inds]
        #TODO maybe: make sure that all voxels of this label correspond to at least one vertex.
        #Create a similar mask for faces by picking only triangles
        #of which all 3 vertices are included
        face_out_mask=np.c_[verts_out_mask[faces[:,0]],verts_out_mask[faces[:,1]],verts_out_mask[faces[:,2]]].all(axis=1)
        #These are the faces to keep...
        faces_out=faces[face_out_mask]
        #...but the old vertices' indexes of faces have to be transformed to the new vrtx_out_inds:
        for iF in range(faces_out.shape[0]):
            for iV in range(3):
                faces_out[iF,iV],= np.where(faces_out[iF,iV]==verts_out_inds)
        #Write the output surfaces to a file
        fsio.write_geometry(out_surf_path,verts_out,faces_out,create_stamp=None,volume_info=volume_info)
        #Create and write output annotations to files
        lab_out=lab[verts_out_inds]
        fs.write_annot(out_surf_path+".annot", lab_out, ctab, names)
        #Write files with the indexes of vertices to keep
        np.save(out_surf_path+"-idx.npy",verts_out_inds)
        np.savetxt(out_surf_path+"-idx.txt",verts_out_inds,fmt='%d')    
            
            
#------------------Subparcellation-subsegmentation-----------------------------    
def make_subparc(v, f, annot, roi_names, ctab, trg_area=100.0):
    # TODO subcort subparc with geodesic on bounding gmwmi
    # TODO normalize fiber counts by relevant gmwmi area

    # build vertex -> face list map
    vfm = [set([]) for _ in v]
    for i, face in enumerate(f):
        for j in face:
            vfm[j].add(i)
    vfm = np.array(vfm)

    # make new annotation
    new_annot = annot.copy()
    new_names = [] # such that new_names[new_annot[i]] is correct name for i'th vertex
    next_aval = 1
    for i in np.unique(annot):
        name = roi_names[i]
        mask = annot == i

        # "unknown", just skip
        if i == -1:
            new_annot[mask] = 0
            new_names.append(name)
            continue

        # indices of faces in ROI
        rfi = set([])
        for face_set in vfm[mask]:
            rfi.update(face_set)
        rfi = np.array(list(rfi))

        # empty roi
        if rfi.size == 0:
            continue

        # compute area of faces in roi
        roi_area = np.sum(tri_area(v[f[rfi]]))

        # choose k for desired roi area
        k = int(roi_area / trg_area) + 1
        assert k>=1

        # cluster centered vertices
        v_roi = v[mask]
        _, i_lab = scipy.cluster.vq.kmeans2(v_roi - v_roi.mean(axis=0), k)

        # update annot
        new_annot[mask] = next_aval + i_lab
        next_aval += k
        new_names += ['%s-%d' % (name.decode('ascii'), j) for j in range(k)]

    # create random colored ctab
    new_ctab = np.random.randint(255, size=(len(new_names), 5))
    r, g, b, _, _ = new_ctab.T
    new_ctab[:, 3] = 0
    new_ctab[:, 4] = r + 256 * g + 256 * 256 * b # fs magic values

    return new_annot, new_ctab, new_names


def subparc_files(hemi, parc_name, out_parc_name, trg_area):
    trg_area = float(trg_area)
    v, f = read_surf(hemi, 'sphere')
    lab, ctab, names = read_annot(hemi, parc_name)
    new_lab, new_ctab, new_names = make_subparc(v, f, lab, names, ctab, trg_area=trg_area)
    write_annot(hemi, out_parc_name, new_lab, new_ctab, new_names)


            
      #d2t=None,       
def connectivity_geodesic_subparc(surf_path,annot_path,con_verts_idx,out_annot_path=None,
                                  ref_vol_path=None,consim_path=None,parc_area=100,
                                  labels=None,hemi=None, mode="con+geod+adj", vox2rastkr_path=None,
                                  lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt')):                                      
    from scipy.spatial.distance import cdist
    from sklearn.cluster import AgglomerativeClustering  
    from scipy.sparse.csgraph import connected_components
    if "geod" in mode:
        from scipy.sparse.csgraph import shortest_path
    #Read the surface...
    (verts, faces,volume_info) = fsio.read_geometry(surf_path,read_metadata=True)
    #...its annotation
    lab, ctab, names = fs.read_annot(annot_path)
    #...and get the correspoding labels:
    labels_annot=annot_names_to_labels(names,hemi,lut_path=lut_path)
    #Read the indexes of vertices neighboring tract ends voxels:
    con_verts_idx=np.load(con_verts_idx)
    #Set the target labels:
    (labels,nLbl)=read_input_labels(labels=labels,hemi=hemi)
    if "con" in mode:  
        #Load voxel connectivity similarity matrix:
        con=np.load(consim_path)
        #Read the DTI to T1 transform: 
        #d2t=np.loadtxt(d2t)
        #Read the reference tdi_lbl volume:
        vollbl=nbl.load(ref_vol_path)
        vox=vollbl.get_data().astype('i')
        voxdim=vox.shape
        #Get only the voxels that correspond to connectome nodes: 
        voxijk,=np.where(vox.flatten()>0)
        voxijk=np.unravel_index(voxijk, voxdim)
        vox=vox[voxijk[0],voxijk[1],voxijk[2]]
        #...and their coordinates in surface ras xyz space
        vox2rastkr=np.loadtxt(vox2rastkr_path)
        voxxzy=vox2rastkr.dot(np.c_[voxijk[0],voxijk[1],voxijk[2], np.ones(vox.shape[0])].T)[:3].T
        #...and transform them to the T1 RAS space of the surface:
        #voxxzy=d2t.dot(np.c_[voxxzy[:,0],voxxzy[:,1],voxxzy[:,2], np.ones(voxxzy.shape[0])].T)[:3].T
        del vollbl, voxijk
    #Initialize the output:
    out_names=[]
    out_ctab=[]
    out_lab=np.array(lab)
    nL=0
    #For every annotation name:
    for iL in range(len(names)):
        print str(iL)+". "+names[iL]
        #Get the mask of the respective vertices
        iVmask=lab==iL
        #...and their indexes
        iV,=np.where(iVmask)
        nV=len(iV)
        #TODO: reconsider this point
        if nV==0:
            continue
        #Find the corresponding label:
        lbl=labels_annot[iL]
        #If it is not one of the target labels:
        if lbl not in labels:
            #Just add this label to the new annotation as it stands:
            out_names.append(names[iL])
            out_ctab.append(ctab[iL])
            #and change the output indices
            nL+=1
            out_lab[iV]=nL
            continue   
        #Get the vertices and faces of this label:
        (verts_lbl,faces_lbl)=extract_subsurf(verts,faces,iVmask)
        #Compute distances among directly connected vertices
        dist=vertex_connectivity(verts_lbl,faces_lbl,mode="sparse",metric='euclidean')
        if "con" in mode:  
            #Mask of label vertices that are neighbors of tract end voxels ("con"):
            iV_con=np.in1d(iV,con_verts_idx)
            nVcon=np.sum(iV_con)
            #Get only the connectome neighboring vertices and faces of this label:
            (verts_lbl_con,faces_lbl_con)=extract_subsurf(verts_lbl,faces_lbl,iV_con)
        else:
            #We work with iV_con from now on even if it is identical to iV
            nVcon=len(iV)
            iV_con=np.ones(nVcon,dtype='bool') 
            verts_lbl_con=verts_lbl
            faces_lbl_con=faces_lbl
        #Calculate total area:
        roi_area=np.sum(tri_area(verts_lbl_con[faces_lbl_con]))
        print "Total ROI area = "+str(roi_area)+" mm2"
        #Calculate number of parcels:
        nParcs = int(np.round(roi_area/parc_area))
        print "Number of parcels for clustering: nParcs="+str(nParcs)
        #If no further parcellation is needed
        if nParcs<2:
            #Just add this label to the new annotation as it stands:
            out_names.append(names[iL])
            out_ctab.append(ctab[iL])
            #and change the output indices
            out_lab[iV]=nL+1
            nL+=1
            continue
        print "Compute geodesic distance matrix"
        geodist=shortest_path(dist.todense(), method='auto', directed=False, return_predecessors=False, unweighted=False, overwrite=False)
        print "Compute structural connectivity constraint matrix"
        connectivity=dist[iV_con,:][:,iV_con]
        connectivity.data=np.ones(connectivity.data.shape)    
        #Check how many connected components there are, and remove all of those with less than 1% of the vertices,
        #by moving them to the noncon group
        (n_components,concomp_labels)=connected_components(connectivity, directed=False, connection='weak', return_labels=True)
        h,_=np.histogram(concomp_labels,np.array(range(n_components+1))-0.5)        
        print 'Before correction: '+str(n_components)+' connected components with '+str(h)+' vertices each'
        #Correction: removal of too small connected components
        #Too small = less than 10% of the expeected average parcel size
        nMinVertsPerComp = int(np.round(0.1*nVcon/nParcs))
        if np.any(h<nMinVertsPerComp):
            for iC in range(n_components):  
                if h[iC]<nMinVertsPerComp:
                    iV_con[concomp_labels==iC] =False
            nVcon=np.sum(iV_con)       
            #Update connectivity matrix after correction if needed:    
            #if "adj" in mode:
            connectivity=dist[iV_con,:][:,iV_con]
            (n_components,concomp_labels)=connected_components(connectivity, directed=False, connection='weak', return_labels=True)
            h,_=np.histogram(concomp_labels,np.array(range(n_components+1))-0.5)        
            print 'After correction: '+str(n_components)+' connected components with '+str(h)+' vertices each'
            #else:
        if "adj" not in mode:    
            connectivity=None    
        del dist    
        if "con" in mode: 
            #Get again after correction the connectome neighboring vertices and faces of this label:
            (verts_lbl_con,faces_lbl_con)=extract_subsurf(verts_lbl,faces_lbl,iV_con)
            print "Compute connectivity similarity affinity matrix..."
            #TODO?: to use aparc+aseg to correspond vertices only to voxels of the same label
            #There would have to be a vertex->voxel of aparc+aseg of the same label -> voxel of tdi_lbl_in_T1 mapping
            #Maybe redundant  because we might be ending to the same voxel of tdi_lbl anyway...
            #Something to test/discuss...
            #Find for each vertex the closest voxel node xyz coordinates:
            v2n = np.argmin(cdist(verts_lbl_con, voxxzy, 'euclidean'),axis=1)
            v2n=vox[v2n]
            print "...whereby vertices correspond to "+str(np.size(np.unique(v2n)))+" distinct voxel nodes"
            #... convert connectivity similarity to a distance matrix
            affinity=1-con[v2n-1,:][:,v2n-1]
            del v2n
        else:
            #Initialize affinity matrix with zeros     
            affinity=np.zeros((nVcon,nVcon)) 
        #Treat the indexes of non-"con" vertices:    
        iV_noncon=~iV_con
        iV_con,=np.where(iV_con)
        if np.any(iV_noncon):
            iV_noncon,=np.where(iV_noncon)
            nVnoncon=len(iV_noncon)
            #Find the closest neighbors of each non-con-vertex to a con-vertex...
            noncon2con_dists=geodist[iV_noncon,:][:,iV_con]
            noncon2con=np.argmin(noncon2con_dists,axis=1)
            noncon2con_distsmins=noncon2con_dists[range(nVnoncon),noncon2con]
            #For every infinite geodesic distance, find the con vertex of minimum euclidean distance
            for iC in range(nVnoncon):
                if np.isinf(noncon2con_distsmins[iC]):
                    noncon2con[iC]=np.argmin(cdist(np.expand_dims(verts_lbl[iV_noncon[iC],:], 1).T, verts_lbl[iV_con,:], 'euclidean'))
            #...and map them  
            noncon2con=iV[noncon2con]
        if "geod" in mode: 
            print "Computing geodesic similarity affinity matrix"
            geodist=geodist[iV_con,:][:,iV_con]
            #Find the maximum non infite geodesic distance:
            temp_ind=np.isfinite(geodist)          
            max_gdist=np.max(geodist[temp_ind],axis=None)
            #Set inf values to maximum distance
            geodist[~temp_ind]=max_gdist
            #Convert them to normalized distances
            geodist=geodist/max_gdist
            #...and add them to the affinity metric
            affinity+=geodist
        del geodist, temp_ind
        print "Run clustering" 
        model = AgglomerativeClustering(n_clusters=nParcs, affinity="precomputed", connectivity=connectivity, linkage='average')
        model.fit(affinity) 
        h,_=np.histogram(model.labels_,np.array(range(nParcs+1))-0.5)        
        print str(nParcs)+' parcels with '+str(h)+' "connectome" vertices each'
        print "Generate new annotations"
        #Create the new annotation for these sub-labels:
        resLbl=[names[iL]+"-"+str(item).zfill(2) for item in np.unique(model.labels_) ]  
        #Add the new label names
        out_names.append(resLbl)
        #Initialize ctab for these labels
        ctab_lbl=np.repeat(ctab[iL,np.newaxis],nParcs,axis=0)
        #For the RGB coordinate with the bigger distance to 255 or 0
        #distribute that distance  to nParcs values:
        iC=np.argsort(ctab[iL,:3])
        x=255-ctab[iL,iC[0]]>=ctab[iL,iC[2]]
        dist=np.where(x,255-ctab[iL,iC[0]],-ctab[iL,iC[2]])
        iC=np.where(x,iC[0],iC[2])
        step=dist/(nParcs-1)
        dist=step*nParcs
        ctab_lbl[:,iC]=np.array(range(ctab[iL,iC],ctab[iL,iC]+dist,step),dtype='int')
        ctab_lbl[:,:3][ctab_lbl[:,:3]<0]=0
        ctab_lbl[:,:3][ctab_lbl[:,:3]>255]=255
        ctab_lbl[:,4]=np.array([rgb_to_fs_magic_number(ctab_lbl[iP,:3]) for iP in range(nParcs)])
        out_ctab.append(ctab_lbl)
        #Finally change the output indices
        #of the "con" vertices...
        out_lab[iV[iV_con]]=nL+model.labels_
        #and of the "non-con" vertices that follow their nearest con neighbor
        if np.any(iV_noncon):
            out_lab[iV[iV_noncon]]=out_lab[noncon2con]
        temp_lbls=np.unique(out_lab[iV[iV_con]])    
        h,hc=np.histogram(out_lab[iV[iV_con]],np.array(np.r_[temp_lbls,temp_lbls[-1]+1])-0.5)  
        for iP in range(nParcs):
            print 'parcel '+str(int(hc[iP]+0.5))+'with '+str(h[iP])+' total vertices'    
        nL+=nParcs
    print "Write output annotation file"
    #Stack everything together
    out_ctab=np.vstack(out_ctab)                        
    out_names=np.hstack(out_names)        
    #...and write the annotation to a file
    if out_annot_path is None:
        out_annot_path=os.path.splitext(annot_path)[0]+str(parc_area)+".annot"
    print out_annot_path
    fs.write_annot(out_annot_path,out_lab, out_ctab, out_names)    


        

         
      

#---------------------------------Dipoles--------------------------------------

def gen_dipole_triplets(pos):
    pos3 = np.repeat(pos, 3, axis=0)
    ori3 = np.tile(np.eye(3), (len(pos), 1))
    return pos3, ori3

def gen_dipoles(pos, ori_or_face=None, out_fname=None):
    "Generate dipoles (or equiv. file) for OpenMEEG."
    if ori_or_face is None:
        pos, ori = gen_dipole_triplets(pos)
    else:
        if ori_or_face.dtype in np.floattypes:
            ori = ori_or_face
        else:
            ori = gen_vertex_normals(v=pos, f=ori_or_face)
    np.savetxt(out_fname, np.c_[pos, ori], fmt='%f')



#-------------------------------Contacts---------------------------------------
def periodic_xyz_for_object(lab, val, aff, bw=0.1, doplot=False):
    "Find blob centers for object in lab volume having value val."
    # TODO handle oblique with multiple spacing
    # vox coords onto first mode
    vox_idx = np.argwhere(lab == val)
    xyz = aff.dot(np.c_[vox_idx, np.ones(vox_idx.shape[0])].T)[:3].T
    xyz_mean = xyz.mean(axis=0)
    xyz -= xyz_mean
    u, s, vt = np.linalg.svd(xyz, 0)
    xi = u[:, 0] * s[0]
    # histogram and ft to find spacing and offset
    bn, bxi_ = np.histogram(xi, np.r_[min(xi) - 0.5 : max(xi) + 0.5 : bw])
    bxi = bxi_[:-1] + bw / 2.0
    w = np.r_[1.0 : 6.0 : 1000j]
    f = (1.0 / w)[:, None]
    Bf = (np.exp(-2 * np.pi * 1j * bxi * f) * bn * bw).sum(axis=-1)
    i_peak = np.argmax(np.abs(Bf))
    theta = np.angle(Bf[i_peak])
    print("[periodic_xyz_for_object]", val, 1/f[i_peak][0], theta)
    xi_o = -theta / (2 * np.pi * f[i_peak])
    xi_pos = np.r_[xi_o : xi.max() : w[i_peak]]
    xi_neg = np.r_[-xi_o : -xi.min() : w[i_peak]]
    xi_pos = np.sort(np.r_[-xi_neg, xi_pos[1:]])
    xyz_pos = np.c_[xi_pos, np.zeros((len(xi_pos), 2))].dot(vt) + xyz_mean
    if doplot:
        pl.figure()
        pl.subplot(2, 1, 1)
        pl.plot(bxi, bn)
        pl.subplot(2, 1, 2)
        pl.plot(w, np.abs(Bf))
        pl.subplot(2, 1, 1)
        cos_arg = 2*np.pi*f[i_peak]*bxi + theta
        pl.plot(bxi, np.cos(cos_arg)*bn.std()+bn.mean(), 'k--', alpha=0.5)
        [pl.axvline(xp, color='r') for xp in xi_pos];
        pl.show()
    return xyz_pos



    

    



if __name__ == '__main__':
    cmd = sys.argv[1]

    if cmd == 'gdist':
        compute_gdist_mat(*sys.argv[2:])
    if cmd == 'subparc':
        subparc_files(*sys.argv[2:])

