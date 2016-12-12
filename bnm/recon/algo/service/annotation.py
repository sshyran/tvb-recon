# -*- coding: utf-8 -*-

import os
import numpy
from collections import OrderedDict
from nibabel.freesurfer import read_annot, write_annot
from bnm.recon.qc.parser.annotation import AnnotationParser


class AnnotationService(object):
    def __init__(self):
        self.annotationParser = AnnotationParser()

    #TODO do we need this?
    def annot_path(self, hemi, annot_name):
        annot_fname = '%s.%s.annot' % (hemi, annot_name)
        annot_path = os.path.join(os.environ['SUBJECTS_DIR'], os.environ['SUBJECT'], 'label', annot_fname)
        return annot_path

    #TODO use this from qc.parser?
    def read_annot(self, hemi, annot_name):
        return read_annot(self.annot_path(hemi, annot_name))

    def write_annot(self, hemi, annot_name, labels, ctab, names):
        return write_annot(self.annot_path(hemi, annot_name), labels, ctab, names)

    def read_lut(self, lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt'), key_mode='label'):
        f = open(lut_path, "r")
        l = list(f)
        f.close()
        ii = -1
        if key_mode == 'label':
            names = OrderedDict()
            colors = OrderedDict()
            labels = []
            for line in l:
                temp = line.split()
                try:
                    label = int(temp[0])
                    ii += 1
                    labels.append(label)
                    names[labels[ii]] = temp[1]
                    colors[labels[ii]] = [int(temp[2]), int(temp[3]), int(temp[4]), int(temp[5])]
                except:
                    pass
        elif key_mode == 'name':
            labels = OrderedDict()
            colors = OrderedDict()
            names = []
            for line in l:
                temp = line.split()
                try:
                    label = int(temp[0])
                    ii += 1
                    names.append(temp[1])
                    labels[names[ii]] = label
                    colors[names[ii]] = [int(temp[2]), int(temp[3]), int(temp[4]), int(temp[5])]
                except:
                    pass
        return (labels, names, colors)

    # TODO where to put this
    def rgb_to_fs_magic_number(self,rgb):
        return rgb[0] + 256 * rgb[1] + 256 * 256 * rgb[2]

    def annot_to_lut(self, hemi, annot_name, lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):
        _, ctab, names = self.read_annot(hemi, annot_name)
        with open(lut_path, 'w') as fd:
            for name, (r, g, b, a, id) in zip(names, ctab):
                fd.write('%d\t%s\t%d %d %d %d\n' % (id, name, r, g, b, a))

    def lut_to_annot_names_ctab(self, lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt'), labels=None):
        (_, names_dict, colors) = self.read_lut(lut_path=lut_path)
        if labels is None:
            labels = names_dict.keys()
        elif isinstance(labels, basestring):
            labels = numpy.array(labels.split()).astype('i').tolist()
        else:
            labels = numpy.array(labels).astype('i').tolist()
        names = []
        ctab = []
        for lbl in labels:
            names.append(names_dict[lbl])
            rgb = numpy.array(colors[lbl])[:3].astype('int64')
            magic_number = self.rgb_to_fs_magic_number(rgb) * numpy.ones((1,), dtype='int64')
            ctab.append(numpy.concatenate([rgb, numpy.zeros((1,), dtype='int64'), magic_number]))
        ctab = numpy.asarray(ctab).astype('int64')
        return (names, ctab)

    def annot_names_to_labels(self, names, ctx=None, lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):
        (labels_dict, _, _) = self.read_lut(lut_path=lut_path, key_mode='name')
        labels = []
        if ctx == 'lh' or ctx == 'rh':
            ctx = 'ctx-' + ctx + '-'
        else:
            ctx = ''
        for name in names:
            labels.append(labels_dict[ctx + name])
        return labels

    def annot_to_conn_conf(self, hemi, annot_name, conn_conf_path):
        _, _, names = self.read_annot(hemi, annot_name)
        with open(conn_conf_path, 'w') as fd:
            for id, name in enumerate(names):
                fd.write('%d\t%s\n' % (id, name))

    def read_input_labels(self, labels=None, hemi=None):
        if labels is not None:
            if isinstance(labels, basestring):
                # Set the target labels
                labels = numpy.array(labels.split()).astype('i').tolist()
        else:
            labels = []
        if hemi is not None:
            hemi = hemi.split()
            for h in hemi:
                if h == 'lh':
                    labels = labels + range(1000, 1036)
                elif h == 'rh':
                    labels = labels + range(2000, 2036)
        nLbl = len(labels)
        return (labels, nLbl)
