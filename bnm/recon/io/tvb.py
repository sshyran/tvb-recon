"""
I/O for TVB formats.

"""

import zipfile
from .core import np_save_strio, StringIO


class TVBWriter(object):
    def write_surface_zip(zip_fname, v, f):
        "Write a surface to a TVB-format ZIP file."
        sv = np_save_strio(v, '%f')
        sf = np_save_strio(f, '%d')
        szf = StringIO()
        zf = zipfile.ZipFile(szf, 'w')
        zf.writestr('vertices.txt', sv.getvalue())
        zf.writestr('triangles.txt', sf.getvalue())
        zf.close()
        with open(zip_fname, 'wb') as fd:
            fd.write(szf.getvalue())
