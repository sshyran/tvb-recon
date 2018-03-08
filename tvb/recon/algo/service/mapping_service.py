# -*- coding: utf-8 -*-
import numpy
from tvb.recon.dax import AtlasSuffix
from tvb.recon.logger import get_logger
from tvb.recon.model.annotation import Annotation


class MappingService(object):
    CORT_TYPE = "aparc"
    SUBCORT_TYPE = "aseg"

    FS_PREFIX_LH_DEFAULT = "ctx-lh-"
    FS_PREFIX_RH_DEFAULT = "ctx-rh-"
    FS_PREFIX_LH_A2009S = "ctx_lh_"
    FS_PREFIX_RH_A2009S = "ctx_rh_"

    UNKNOWN_REGION = "unknown"
    CORPUSCALLOSUM_REGION = "corpuscallosum"
    UNKNOWN_SUBCORTICAL_REGION = "Unknown"
    CORPUSCALLOSUM_REGION_A2009S = "Corpus_callosum"

    atlas_suffix = AtlasSuffix.DEFAULT
    fs_prefix_lh = FS_PREFIX_LH_DEFAULT
    fs_prefix_rh = FS_PREFIX_RH_DEFAULT

    logger = get_logger(__name__)

    def __init__(self, atlas_suffix: AtlasSuffix, cort_annot_lh: Annotation, cort_annot_rh: Annotation,
                 subcort_annot_lh: Annotation, subcort_annot_rh: Annotation):
        self.atlas_suffix = atlas_suffix
        if atlas_suffix == AtlasSuffix.A2009S:
            self.fs_prefix_lh = self.FS_PREFIX_LH_A2009S
            self.fs_prefix_rh = self.FS_PREFIX_RH_A2009S

        self.cort_lut_dict = self.generate_lut_dict_from_annot(cort_annot_lh, cort_annot_rh, self.CORT_TYPE, 0)
        self.subcort_lut_dict = self.generate_lut_dict_from_annot(subcort_annot_lh, subcort_annot_rh,
                                                                  self.SUBCORT_TYPE, len(self.cort_lut_dict))
        self.cort_region_mapping = list()
        self.subcort_region_mapping = list()

    def get_lh_regions(self):
        return list(self.cort_lut_dict.keys())[0:self.len_lh]

    def get_rh_regions(self):
        return list(self.cort_lut_dict.keys())[self.len_rh:]

    def generate_lut_dict_from_annot(self, annot_lh: Annotation, annot_rh: Annotation, annot_type: str,
                                     idx: int) -> dict:
        dict_lh = self._get_dict_from_annot(annot_lh)
        dict_rh = self._get_dict_from_annot(annot_rh)

        if annot_type == self.CORT_TYPE:
            self.len_lh = len(dict_lh)
            self.len_rh = len(dict_rh)
            return self._prepare_cort_lut_dict(dict_lh, dict_rh, idx)

        return self._prepare_subcort_lut_dict(dict_lh, dict_rh, idx)

    def _get_dict_from_annot(self, annot: Annotation) -> dict:
        annot_dict = dict()
        vtx_rm = annot.region_mapping
        vtx_rm_unique_vals = numpy.unique(vtx_rm)

        region_names = annot.region_names
        for idx, region_name in enumerate(region_names):
            if self.atlas_suffix == AtlasSuffix.A2009S:
                if "&" in region_name:
                    region_names[idx] = region_name.replace("&", "_and_")

        for unwanted_region in (self.UNKNOWN_REGION, self.UNKNOWN_SUBCORTICAL_REGION, self.CORPUSCALLOSUM_REGION,
                                self.CORPUSCALLOSUM_REGION_A2009S):
            if unwanted_region in region_names and region_names.index(unwanted_region) in vtx_rm_unique_vals:
                self.logger.warn("This annotation contains vertices for %s" % unwanted_region)

        region_names_to_keep = [region_names[idx] for idx in range(len(region_names)) if idx in vtx_rm_unique_vals]

        outside_range_values = list(set(vtx_rm_unique_vals) - set(range(len(region_names))))
        if len(outside_range_values) > 0:
            self.logger.warn("This annotation contains vertices associated to values outside the interval [ %d, %d]"
                             % (0, len(region_names) - 1))
            self.logger.warn("These values are: %s" % outside_range_values)
            self.logger.info("Vertices mapped to these values will be mapped to %s" % self.UNKNOWN_REGION)
            if self.atlas_suffix == AtlasSuffix.A2009S:
                region_names_to_keep.append(self.UNKNOWN_SUBCORTICAL_REGION)
            else:
                region_names_to_keep.append(self.UNKNOWN_REGION)

        for idx, name in enumerate(region_names_to_keep):
            annot_dict[idx] = name

        return annot_dict

    def _prepare_cort_lut_dict(self, dict_lh: dict, dict_rh: dict, idx: int) -> dict:
        lut_dict = dict()

        lut_dict.update({idx + key: self.fs_prefix_lh + val for (key, val) in dict_lh.items()})
        idx += len(lut_dict)
        lut_dict.update({idx + key: self.fs_prefix_rh + val for (key, val) in dict_rh.items()})

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

        lh_annot.region_mapping[lh_annot.region_mapping == -1] = 0
        rh_annot.region_mapping[rh_annot.region_mapping == -1] = 0

        self.lh_region_mapping = list()
        for lbl in lh_annot.region_mapping:
            current_region_name = lh_annot.region_names[lbl]
            region_mapping.append(cort_inv_lut_dict.get(self.fs_prefix_lh + current_region_name))
            self.lh_region_mapping.append(cort_inv_lut_dict.get(self.fs_prefix_lh + current_region_name))

        self.rh_region_mapping = list()
        for lbl in rh_annot.region_mapping:
            current_region_name = rh_annot.region_names[lbl]
            region_mapping.append(cort_inv_lut_dict.get(self.fs_prefix_rh + current_region_name))
            self.rh_region_mapping.append(cort_inv_lut_dict.get(self.fs_prefix_rh + current_region_name))

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
            src_ind = lut_idx_to_name_dict.get(trg_name, None)
            if src_ind is not None:
                if trg_name == self.UNKNOWN_SUBCORTICAL_REGION:
                    self.logger.warn(
                        "The subcortical surfaces contain the region: %s and it will be mapped to -1 for aparc+aseg" % self.UNKNOWN_SUBCORTICAL_REGION)
                    src_to_trg[src_ind] = -1
                else:
                    src_to_trg[src_ind] = trg_ind

        return src_to_trg

    def get_mapping_for_connectome_generation(self):
        non_zero_keys_dict = ({key + 1: val for (key, val) in self.get_entire_lut().items()})
        return non_zero_keys_dict
