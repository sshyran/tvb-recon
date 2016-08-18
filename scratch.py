import pyqtgraph as pg
import pyqtgraph.opengl as gl
import nibabel.freesurfer
import scipy.cluster


def plot_points(v, size=1.5, title='', labels=None):
    win = gl.GLViewWidget()
    win.show()
    win.setWindowTitle(title)
    item = gl.GLScatterPlotItem(pos=v, size=size)
    if labels is not None:
        rc = np.random.rand(np.unique(labels).size, 4) / 2.0 + 0.5
        item.setData(color=rc[labels])
    win.addItem(item)
    return win, item


def plot_label(idx, v, v_labels):
    return plot_points(v[idx == v_labels])


def sublab_roi(idx, v, v_labels, k):
    "Cluster vertices in given label into k clusters."
    v_roi = v[idx == v_labels]
    v_roi -= v_roi.mean(axis=0)
    _, sublab = scipy.cluster.vq.kmeans2(v_roi, k)
    return v_roi, sublab


def plot_sublabel(idx, v, v_labels, k):
    v_roi, sublab = sublab_roi(idx, v, v_labels, k)
    return plot_points(v_roi, labels=sublab)


def vert_face_map(v, f):
    "Compute mapping from vertix index to set of face indices."
    faces_for_vertex = [set([]) for _ in v]
    for i, face in enumerate(f):
        for j in face:
            faces_for_vertex[j].add(i)
    return np.array(faces_for_vertex)


def roi_faces(idx, vfm, v_labels):
    "Generate aoi. of faces w/ 1+ vertex in label."
    total_face_set = set([])
    for face_set in vfm[v_labels == idx]:
        total_face_set.update(face_set)
    return np.array(list(total_face_set))


def roi_area(idx, vfm, v_labels):
    rfi = roi_faces(idx, vfm, v_labels)
    if rfi.size == 0:
        return 0.0
    tri_xyz = v[f[rfi]]
    i, j, k = transpose(tri_xyz, (1, 0, 2))
    ij = j - i
    ik = k - i
    return np.sum(np.sqrt(np.sum(np.cross(ij, ik)**2, axis=1)) / 2.0)


path = '/usr/local/freesurfer/subjects/junch/'
v, f = nibabel.freesurfer.read_geometry(path + 'surf/lh.white')
v_roi, ctab, roi_names = nibabel.freesurfer.read_annot(path + 'label/lh.aparc.annot')
