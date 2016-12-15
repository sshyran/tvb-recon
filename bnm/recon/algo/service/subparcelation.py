# -*- coding: utf-8 -*-

import os
import numpy
import nibabel
import scipy
from scipy.spatial.distance import cdist, squareform
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster import hierarchy
from scipy.sparse.csgraph import connected_components
from scipy.sparse.csgraph import shortest_path
from nibabel.freesurfer.io import read_geometry, read_annot, write_annot
from bnm.recon.algo.service.annotation import AnnotationService
from bnm.recon.algo.service.surface import SurfaceService

class SubparcellationService(object):
    def __init__(self):
        self.annotationService = AnnotationService()
        self.surfaceService = SurfaceService()

    def make_subparc(self, v, f, annot, roi_names, ctab, trg_area=100.0):
        # TODO subcort subparc with geodesic on bounding gmwmi
        # TODO normalize fiber counts by relevant gmwmi area

        # build vertex -> face list map
        vfm = [set([]) for _ in v]
        for i, face in enumerate(f):
            for j in face:
                vfm[j].add(i)
        vfm = numpy.array(vfm)

        # make new annotation
        new_annot = annot.copy()
        new_names = []  # such that new_names[new_annot[i]] is correct name for i'th vertex
        next_aval = 1
        for i in numpy.unique(annot):
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
            rfi = numpy.array(list(rfi))

            # empty roi
            if rfi.size == 0:
                continue

            # compute area of faces in roi
            roi_area = numpy.sum(self.surfaceService.tri_area(v[f[rfi]]))

            # choose k for desired roi area
            k = int(roi_area / trg_area) + 1
            assert k >= 1

            # cluster centered vertices
            v_roi = v[mask]
            _, i_lab = scipy.cluster.vq.kmeans2(v_roi - v_roi.mean(axis=0), k)

            # update annot
            new_annot[mask] = next_aval + i_lab
            next_aval += k
            new_names += ['%s-%d' % (name.decode('ascii'), j) for j in range(k)]

        # create random colored ctab
        new_ctab = numpy.random.randint(255, size=(len(new_names), 5))
        r, g, b, _, _ = new_ctab.T
        new_ctab[:, 3] = 0
        new_ctab[:, 4] = r + 256 * g + 256 * 256 * b  # fs magic values

        return new_annot, new_ctab, new_names

    def subparc_files(self, hemi, parc_name, out_parc_name, trg_area):
        trg_area = float(trg_area)
        v, f = self.surfaceService.read_surf(hemi, 'sphere')
        lab, ctab, names = read_annot(hemi, parc_name)
        new_lab, new_ctab, new_names = self.make_subparc(v, f, lab, names, ctab, trg_area=trg_area)
        write_annot(self.annotationService.annot_path(hemi, out_parc_name), new_lab, new_ctab, new_names)

        # d2t=None,
    def connectivity_geodesic_subparc(self, surf_path, annot_path, con_verts_idx, out_annot_path=None,
                                      parc_area=100, labels=None, hemi=None, ctx=None, mode="con+geod+adj",
                                      cras_path=None, ref_vol_path=None, consim_path=None,
                                      lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):

        # Read the surface...
        (verts, faces, volume_info) = read_geometry(surf_path, read_metadata=True)
        # ...its annotation
        lab, ctab, names = read_annot(annot_path)
        # ...and get the correspoding labels:
        labels_annot = self.annotationService.annot_names_to_labels(names, ctx, lut_path=lut_path)
        # Read the indexes of vertices neighboring tract ends voxels:
        con_verts_idx = numpy.load(con_verts_idx)
        # Set the target labels:
        (labels, nLbl) = self.annotationService.read_input_labels(labels=labels, hemi=hemi)
        if "con" in mode:
            # Load voxel connectivity similarity matrix:
            con = numpy.load(consim_path).astype('single')
            # Read the cras:
            cras = numpy.loadtxt(cras_path)
            # Read the DTI to T1 transform:
            # d2t=numpy.loadtxt(d2t)
            # Read the reference tdi_lbl volume:
            vollbl = nibabel.load(ref_vol_path)
            vox = vollbl.get_data().astype('i')
            voxdim = vox.shape
            # Get only the voxels that correspond to connectome nodes:
            voxijk, = numpy.where(vox.flatten() > 0)
            voxijk = numpy.unravel_index(voxijk, voxdim)
            vox = vox[voxijk[0], voxijk[1], voxijk[2]]
            # ...and their coordinates in surface ras xyz space
            vox2ras = vollbl.affine
            voxxzy = vox2ras.dot(numpy.c_[voxijk[0], voxijk[1], voxijk[2], numpy.ones(vox.shape[0])].T)[:3].T
            # ...and transform them to the T1 RAS space of the surface:
            # voxxzy=d2t.dot(numpy.c_[voxxzy[:,0],voxxzy[:,1],voxxzy[:,2], numpy.ones(voxxzy.shape[0])].T)[:3].T
            del vollbl, voxijk
        # Initialize the output:
        out_names = []
        out_ctab = []
        out_lab = numpy.array(lab)
        nL = -1
        # For every annotation name:
        for iL in range(len(names)):
            print str(iL) + ". " + names[iL]
            # Get the mask of the respective vertices
            iVmask = lab == iL
            # ...and their indexes
            iV, = numpy.where(iVmask)
            nV = len(iV)
            # Find the corresponding label:
            lbl = labels_annot[iL]
            # If there are no associated vertices or if it is not one of the target labels:
            if nV == 0 or (lbl not in labels):
                # Just add this label to the new annotation as it stands:
                out_names.append(names[iL])
                out_ctab.append(ctab[iL])
                # and change the output indices
                nL += 1
                out_lab[iV] = nL
                continue
                # Get the vertices and faces of this label !in tk ras coordinates!:
            (verts_lbl, faces_lbl) = self.surfaceService.extract_subsurf(verts, faces, iVmask)
            # Compute distances among directly connected vertices
            dist = self.surfaceService.vertex_connectivity(verts_lbl, faces_lbl, mode="sparse", metric='euclidean').astype('single')
            if "con" in mode:
                # Mask of label vertices that are neighbors of tract end voxels ("con"):
                iV_con = numpy.in1d(iV, con_verts_idx)
                nVcon = numpy.sum(iV_con)
                # Get only the connectome neighboring vertices and faces of this label:
                (verts_lbl_con, faces_lbl_con) = self.surfaceService.extract_subsurf(verts_lbl, faces_lbl, iV_con)
            else:
                # We work with iV_con from now on even if it is identical to iV
                nVcon = len(iV)
                iV_con = numpy.ones(nVcon, dtype='bool')
                verts_lbl_con = verts_lbl
                faces_lbl_con = faces_lbl
            # Calculate total area:
            roi_area = numpy.sum(self.surfaceService.tri_area(verts_lbl_con[faces_lbl_con]))
            print "Total ROI area = " + str(roi_area) + " mm2"
            # Calculate number of parcels:
            nParcs = int(numpy.round(roi_area / parc_area))
            print "Number of parcels for clustering: nParcs=" + str(nParcs)
            # If no further parcellation is needed
            if nParcs < 2:
                # Just add this label to the new annotation as it stands:
                out_names.append(names[iL])
                out_ctab.append(ctab[iL])
                # and change the output indices
                nL += 1
                out_lab[iV] = nL
                continue
            print "Compute geodesic distance matrix"
            geodist = shortest_path(dist.todense(), method='auto', directed=False, return_predecessors=False,
                                    unweighted=False, overwrite=False).astype('single')
            print "Compute structural connectivity constraint matrix"
            connectivity = dist[iV_con, :][:, iV_con].astype('single')
            connectivity.data = numpy.ones(connectivity.data.shape).astype('single')
            # Check how many connected components there are, and remove all of those with less than 1% of the vertices,
            # by moving them to the noncon group
            (n_components, concomp_labels) = connected_components(connectivity, directed=False, connection='weak',
                                                                  return_labels=True)
            h, _ = numpy.histogram(concomp_labels, numpy.array(range(n_components + 1)) - 0.5)
            print 'Before correction: ' + str(n_components) + ' connected components with ' + str(h) + ' vertices each'
            # Correction: removal of too small connected components
            # Too small = less than 10% of the expeected average parcel size
            nMinVertsPerComp = int(numpy.round(0.1 * nVcon / nParcs))
            if numpy.any(h < nMinVertsPerComp):
                for iC in range(n_components):
                    if h[iC] < nMinVertsPerComp:
                        iV_con[concomp_labels == iC] = False
                nVcon = numpy.sum(iV_con)
                # Update connectivity matrix after correction if needed:
                # if "adj" in mode:
                connectivity = dist[iV_con, :][:, iV_con].astype('single')
                connectivity.data = numpy.ones(connectivity.data.shape).astype('single')
                (n_components, concomp_labels) = connected_components(connectivity, directed=False, connection='weak',
                                                                      return_labels=True)
                h, _ = numpy.histogram(concomp_labels, numpy.array(range(n_components + 1)) - 0.5)
                print 'After correction: ' + str(n_components) + ' connected components with ' + str(
                    h) + ' vertices each'
                # else:
            if "adj" not in mode:
                connectivity = None
            del dist
            if "con" in mode:
                # Get again after correction the connectome neighboring vertices and faces of this label:
                (verts_lbl_con, faces_lbl_con) = self.surfaceService.extract_subsurf(verts_lbl, faces_lbl, iV_con)
                # Add the cras to take them to scanner ras coordinates:
                verts_lbl_con += numpy.repeat(numpy.expand_dims(cras, 1).T, nVcon, axis=0)
                print "Compute connectivity similarity affinity matrix..."
                # TODO?: to use aparc+aseg to correspond vertices only to voxels of the same label
                # There would have to be a vertex->voxel of aparc+aseg of the same label -> voxel of tdi_lbl_in_T1 mapping
                # Maybe redundant  because we might be ending to the same voxel of tdi_lbl anyway...
                # Something to test/discuss...
                # Find for each vertex the closest voxel node xyz coordinates:
                v2n = numpy.argmin(cdist(verts_lbl_con, voxxzy, 'euclidean'), axis=1)
                v2n = vox[v2n]
                print "...whereby vertices correspond to " + str(numpy.size(numpy.unique(v2n))) + " distinct voxel nodes"
                # ... convert connectivity similarity to a distance matrix
                #affinity = 1 - con[v2n - 1, :][:, v2n - 1]
                affinity = numpy.arccos(con[v2n - 1, :][:, v2n - 1])
                #Normalize:
                affinity /=numpy.max(affinity).astype('single')
                del v2n
            else:
                # Initialize affinity matrix with zeros
                affinity = numpy.zeros((nVcon, nVcon)).astype('single')
                # Treat the indexes of non-"con" vertices:
            iV_noncon = ~iV_con
            iV_con, = numpy.where(iV_con)
            if numpy.any(iV_noncon):
                iV_noncon, = numpy.where(iV_noncon)
                nVnoncon = len(iV_noncon)
                # Find the closest neighbors of each non-con-vertex to a con-vertex...
                noncon2con_dists = geodist[iV_noncon, :][:, iV_con]
                noncon2con = numpy.argmin(noncon2con_dists, axis=1)
                noncon2con_distsmins = noncon2con_dists[range(nVnoncon), noncon2con]
                # For every infinite geodesic distance, find the con vertex of minimum euclidean distance
                for iC in range(nVnoncon):
                    if numpy.isinf(noncon2con_distsmins[iC]):
                        noncon2con[iC] = numpy.argmin(
                            cdist(numpy.expand_dims(verts_lbl[iV_noncon[iC], :], 1).T, verts_lbl[iV_con, :], 'euclidean'))
                # ...and map them
                noncon2con = iV[noncon2con]
            if "geod" in mode:
                print "Computing geodesic similarity affinity matrix"
                geodist = geodist[iV_con, :][:, iV_con]
                # Find the maximum non infite geodesic distance:
                max_gdist = numpy.max(geodist, axis=None)
                if ~numpy.isfinite(max_gdist):
                    print "Non-finite geodesic distance detected. Evidence for non-connected surface."
                    return
                # Convert them to normalized distances
                geodist = geodist / max_gdist
                # ...and add them to the affinity metric
                affinity += geodist
            del geodist   
            print "Run clustering"
            algo="scikit" #or "scipy"
            #TODO: probably we either need or own classification algorithm, or
            #an algorithm for going through the tree returned by scikit and, 
            #re-mering/splitting until we creat clusters within some min-max limits of
            #either number of vertices or total surface area...
            if algo == "scikit":
            #scikit learn algorithm:
                model = AgglomerativeClustering(n_clusters=nParcs, affinity="precomputed", 
                                                connectivity=connectivity, linkage='average')
                model.fit(affinity)
                clusters=model.labels_
                	
                 #You can also do 
                 #children to get an (nVcon,2) array of all nodes with their children
                 #or
                 #tree=dict(enumerate(model.children_, model.n_leaves_)), 
                 #which will give you a dictionary where the each key is the 
                 #ID of a node and the value is the pair of IDs of its children.
                 #or to get a list:
                 #import itertools
                 #ii = itertools.count(nVcon)
                 #tree=[{'node_id': next(ii), 'left': x[0], 'right':x[1]} for x in model.children_]
            else:
                #scipy algorithm:
                numpy.fill_diagonal(affinity,0.0)
                affinity=squareform(affinity).astype('single')
                #affinity[affinity==0.0]=0.001
                Z=hierarchy.ward(affinity)
                try:
                    clusters=hierarchy.fcluster(Z,nParcs,criterion='maxclust')-1
                except:
                    print "Shit!"
                    return
            clusters_labels=numpy.unique(clusters)
            h, _ = numpy.histogram(clusters, numpy.array(clusters_labels) - 0.5)
            print str(nParcs) + ' parcels with ' + str(h) + ' "connectome" vertices each'
            print "Generate new annotations"
            # Create the new annotation for these sub-labels:
            resLbl = [names[iL] + "-" + str(item).zfill(2) for item in clusters_labels]
            # Add the new label names
            out_names.append(resLbl)
            # Initialize ctab for these labels
            ctab_lbl = numpy.repeat(ctab[iL, numpy.newaxis], nParcs, axis=0)
            # For the RGB coordinate with the bigger distance to 255 or 0
            # distribute that distance  to nParcs values:
            iC = numpy.argsort(ctab[iL, :3])
            x = 255 - ctab[iL, iC[0]] >= ctab[iL, iC[2]]
            dist = numpy.where(x, 255 - ctab[iL, iC[0]], -ctab[iL, iC[2]])
            iC = numpy.where(x, iC[0], iC[2])
            step = dist / (nParcs - 1)
            dist = step * nParcs
            ctab_lbl[:, iC] = numpy.array(range(ctab[iL, iC], ctab[iL, iC] + dist, step), dtype='int')
            ctab_lbl[:, :3][ctab_lbl[:, :3] < 0] = 0
            ctab_lbl[:, :3][ctab_lbl[:, :3] > 255] = 255
            ctab_lbl[:, 4] = numpy.array([self.annotationService.rgb_to_fs_magic_number(ctab_lbl[iP, :3]) for iP in range(nParcs)])
            out_ctab.append(ctab_lbl)
            # Finally change the output indices
            # of the "con" vertices...
            out_lab[iV[iV_con]] = nL + clusters
            # and of the "non-con" vertices that follow their nearest con neighbor
            if numpy.any(iV_noncon):
                out_lab[iV[iV_noncon]] = out_lab[noncon2con]
            temp_lbls = numpy.unique(out_lab[iV[iV_con]])
            h, hc = numpy.histogram(out_lab[iV[iV_con]], numpy.array(numpy.r_[temp_lbls, temp_lbls[-1] + 1]) - 0.5)
            #+ str(int(hc[iP] + 0.5))             
            for iP in range(nParcs):
                print 'parcel '+resLbl[iP]+ 'with ' + str(h[iP]) + ' total vertices'
            nL += nParcs
        print "Write output annotation file"
        # Stack everything together
        out_ctab = numpy.vstack(out_ctab)
        out_names = numpy.hstack(out_names)
        # ...and write the annotation to a file
        if out_annot_path is None:
            out_annot_path = os.path.splitext(annot_path)[0] + str(parc_area) + ".annot"
        print out_annot_path
        write_annot(out_annot_path, out_lab, out_ctab, out_names)
