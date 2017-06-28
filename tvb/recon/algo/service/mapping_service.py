# -*- coding: utf-8 -*-
import numpy
from tvb.recon.model.annotation import Annotation


class MappingService(object):
    CORT_TYPE = "aparc"
    SUBCORT_TYPE = "aseg"

    fs_prefix_lh = "ctx-lh-"
    fs_prefix_rh = "ctx-rh-"
    unknown_region = "unknown"

    def __init__(self, cort_annot_lh: Annotation, cort_annot_rh: Annotation, subcort_annot_lh: Annotation,
                 subcort_annot_rh: Annotation):
        self.cort_lut_dict = self.generate_lut_dict_from_annot(cort_annot_lh, cort_annot_rh, self.CORT_TYPE, 0)
        self.subcort_lut_dict = self.generate_lut_dict_from_annot(subcort_annot_lh, subcort_annot_rh, self.SUBCORT_TYPE,
                                                                  len(self.cort_lut_dict))
        self.cort_region_mapping = list()
        self.subcort_region_mapping = list()

    def generate_lut_dict_from_annot(self, annot_lh: Annotation, annot_rh: Annotation, annot_type: str,
                                     idx: int) -> dict:
        dict_lh = self._get_dict_from_annot(annot_lh)
        dict_rh = self._get_dict_from_annot(annot_rh)

        if annot_type == self.CORT_TYPE:
            return self._prepare_cort_lut_dict(dict_lh, dict_rh, idx)

        return self._prepare_subcort_lut_dict(dict_lh, dict_rh, idx)

    def _get_dict_from_annot(self, annot: Annotation) -> dict:
        annot_dict = dict()
        for idx, name in enumerate(annot.region_names):
            annot_dict[idx] = name

        return annot_dict

    def _prepare_cort_lut_dict(self, dict_lh: dict, dict_rh: dict, idx: int) -> dict:
        lut_dict = dict()
        has_cc_anterior = False
        lut_dict[0] = self.unknown_region
        for (key, val) in dict_lh.items():
            if (val != self.unknown_region):
                if (val == "corpuscallosum"):
                    has_cc_anterior = True
                    continue
                if has_cc_anterior:
                    lut_dict.update({idx + key - 1: self.fs_prefix_lh + val})
                else:
                    lut_dict.update({idx + key: self.fs_prefix_lh + val})

        has_cc_anterior = False
        idx += len(lut_dict)
        for (key, val) in dict_rh.items():
            if (val != self.unknown_region):
                if (val == "corpuscallosum"):
                    has_cc_anterior = True
                    continue
                if has_cc_anterior:
                    lut_dict.update({idx + key - 2: self.fs_prefix_rh + val})
                else:
                    lut_dict.update({idx + key - 1: self.fs_prefix_rh + val})

        return lut_dict

    def _prepare_subcort_lut_dict(self, dict_lh: dict, dict_rh: dict, idx: int) -> dict:
        lut_dict = dict()
        lut_dict.update({idx + key: val for (key, val) in dict_lh.items()})

        idx += len(lut_dict)
        lut_dict.update({idx + key: val for (key, val) in dict_rh.items()})

        return lut_dict

    def _invert_color_lut(self, color_lut_dict: dict) -> dict:
        inv_dict = {}
        for (key, val) in color_lut_dict.items():
            inv_dict.update({val: key})
        return inv_dict

    def generate_region_mapping_for_cort_annot(self, lh_annot: Annotation, rh_annot: Annotation):
        region_mapping = list()
        cort_inv_lut_dict = self._invert_color_lut(self.cort_lut_dict)

        for lbl in lh_annot.region_mapping:
            current_region_name = lh_annot.region_names[lbl]
            if current_region_name == self.unknown_region:
                region_mapping.append(cort_inv_lut_dict.get(current_region_name))
            else:
                if (current_region_name == "corpuscallosum"):
                    region_mapping.append(0)
                else:
                    region_mapping.append(cort_inv_lut_dict.get(self.fs_prefix_lh + current_region_name))

        for lbl in rh_annot.region_mapping:
            current_region_name = rh_annot.region_names[lbl]
            if current_region_name == self.unknown_region:
                region_mapping.append(cort_inv_lut_dict.get(current_region_name))
            else:
                if (current_region_name == "corpuscallosum"):
                    region_mapping.append(0)
                else:
                    region_mapping.append(cort_inv_lut_dict.get(self.fs_prefix_rh + current_region_name))

        self.cort_region_mapping = region_mapping

    def generate_region_mapping_for_subcort_annot(self, lh_annot: Annotation, rh_annot: Annotation):
        region_mapping = list()
        subcort_inv_lut_dict = self._invert_color_lut(self.subcort_lut_dict)

        for annot in (lh_annot, rh_annot):
            for lbl in annot.region_mapping:
                region_mapping.append(subcort_inv_lut_dict.get(annot.region_names[lbl]))

        self.subcort_region_mapping = region_mapping

    def is_cortical_region_mapping(self):
        return list(numpy.ones(len(self.cort_lut_dict), dtype=int)) + list(numpy.zeros(len(self.subcort_lut_dict),
                                                                                       dtype=int))

    def get_all_regions(self):
        return list(self.cort_lut_dict.keys()) + list(self.subcort_lut_dict.keys())

    def get_entire_lut(self):
        dict = {}
        dict.update(self.cort_lut_dict)
        dict.update(self.subcort_lut_dict)
        return dict

    def get_mapping_for_aparc_aseg(self, lut_idx_to_name_dict: dict) -> dict:
        trg_names_labels_dict = self._invert_color_lut(self.cort_lut_dict)
        trg_names_labels_dict.update(self._invert_color_lut(self.subcort_lut_dict))

        src_to_trg = dict()
        for trg_name, trg_ind in trg_names_labels_dict.items():
            if trg_name is self.unknown_region:
                src_to_trg[0] = -1
            else:
                src_ind = lut_idx_to_name_dict.get(trg_name, None)
                if src_ind is not None:
                    src_to_trg[src_ind] = trg_ind

        return src_to_trg
