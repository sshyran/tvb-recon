from enum import Enum


class Inputs(Enum):
    T1_INPUT = "t1_input.nii.gz"
    DWI_INPUT = "dwi_input.mif"
    FS_LUT = "fs_color_lut.txt"
    FS_DEFAULT = "fs_default.txt"
    T2_INPUT = "t2_input.nii.gz"
    FLAIR_INPUT = "flair_input.nii.gz"


class T1Files(Enum):
    T1_INPUT_CONVERTED = "t1_input.mgz"
    T2_CONVERTED = "t2.nii.gz"
    FLAIR_CONVERTED = "flair.nii.gz"
    T1_MGZ = "T1.mgz"
    APARC_ASEG_MGZ = "aparc+aseg.mgz"
    T1_NII_GZ = "T1.nii.gz"
    APARC_ASEG_NII_GZ = "aparc+aseg.nii.gz"


class T1JobNames(Enum):
    MRI_CONVERT = "mri_convert"
    RECON_ALL = "recon"


class DWIFiles(Enum):
    DWI_RAW_MIF = "dwi_raw.mif"
    DWI_MIF = "dwi.mif"
    MASK_MIF = "mask.mif"
    B0_NII_GZ = "b0.nii.gz"


class DWIJobNames(Enum):
    MRCONVERT = "mrconvert"
    DWIPREPROC = "dwipreproc"
    DWI2MASK = "dwi2mask"
    DWIEXTRACT = "dwiextract"


class CoregFiles(Enum):
    D2T_MAT = "d2t.mat"
    T2D_MAT = "t2d.mat"
    B0_IN_T1 = "b0-in-t1.nii.gz"
    T1_IN_D = "t1-in-d.nii.gz"
    APARC_AGEG_IN_D = "aparc+aseg-in-d.nii.gz"


class CoregJobNames(Enum):
    FLIRT = "flirt"
    CONVERT_XFM = "convert-xfm"
    FLIRT_REVERSED = "flirt-reversed"


class TractsGenFiles(Enum):
    FILE_5TT_MIF = "5tt.mif"
    GMWMI_MIF = "gmwmi.mif"
    FILE_5TTVIS_MIF = "5ttvis.mif"
    RESPONSE_TXT = "response.txt"
    WM_FOD_MIF = "wm_fod.mif"
    FILE_25M_TCK = "25M.tck"
    FILE_5M_TCK = "5M.tck"
    TDI_ENDS_MIF = "tdi_ends.mif"
    GM_MIF = "gm.mif"
    CSF_MIF = "csf.mif"
    RF_WM = "RF_WM.txt"
    RF_GM = "RF_GM.txt"
    RF_CSF = "RF_CSF.txt"
    RF_VOXELS = "RF_voxels.mif"


class TractsGenJobNames(Enum):
    JOB_5TTGEN = "ttgen"
    JOB_5TT2GMWMI = "tt2gmwmi"
    JOB_5TT2VIS = "tt2vis"
    DWI2RESPONSE = "dwi2response"
    DWI2FOD = "dwi2fod"
    TCKGEN = "tckgen"
    TCKSIFT = "tcksift"
    TCKMAP = "tckmap"
    LABEL_CONVERT = "labelconvert"
    TCK2CONNECTOME = "tck2connectome"
    DWI2RESPONSE_MSMT = "dwi2response-msmt"
    MSDWI2FOD = "msdwi2fod"
