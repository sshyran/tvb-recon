
"""
CLI information for mtrix.

"""

from .core import BaseCLI, BaseEnv, BaseFlags


class BaseMtrixFlags(BaseFlags):
    """
    Base flags for mtrix commands.

    """
    pass


class BaseMtrixEnv(BaseEnv):
    """
    Base environment variables for mtrix commands.

    """
    pass


class BaseMtrixCLI(BaseCLI):
    """
    Base CLI for mtrix commands.

    """

    class Flags(BaseMtrixFlags):
        pass

    class Env(BaseMtrixEnv):
        pass



class fttgen(BaseMtrixCLI):
    """
    The 5ttgen command from the mtrix package.

    """
    exe = '5ttgen'


class ftt2gmwmi(BaseMtrixCLI):
    """
    The 5tt2gmwmi command from the mtrix package.

    """
    exe = '5tt2gmwmi'


class ftt2vis(BaseMtrixCLI):
    """
    The 5tt2vis command from the mtrix package.

    """
    exe = '5tt2vis'


class dwi2fod(BaseMtrixCLI):
    """
    The dwi2fod command from the mtrix package.

    """
    exe = 'dwi2fod'


class dwi2mask(BaseMtrixCLI):
    """
    The dwi2mask command from the mtrix package.

    """
    exe = 'dwi2mask'


class dwi2response(BaseMtrixCLI):
    """
    The dwi2response command from the mtrix package.

    """
    exe = 'dwi2response'


class dwipreproc(BaseMtrixCLI):
    """
    The dwipreproc command from the mtrix package.

    """
    exe = 'dwipreproc'


class dwiextract(BaseMtrixCLI):
    """
    The dwiextract command from the mtrix package.

    """
    exe = 'dwiextract'

    class Flags(BaseMtrixCLI.Flags):
        bzero = '-bzero'


class labelconvert(BaseMtrixCLI):
    """
    The labelconvert command from the mtrix package.

    """
    exe = 'labelconvert'


class mrview(BaseMtrixCLI):
    """
    The mrview command from the mtrix package.

    """
    exe = 'mrview'


class mrconvert(BaseMtrixCLI):
    """
    The mrconvert command from the mtrix package.

    """
    exe = 'mrconvert'


class msdwi2fod(BaseMtrixCLI):
    """
    The msdwi2fod command from the mtrix package.

    """
    exe = 'msdwi2fod'


class tck2connectome(BaseMtrixCLI):
    """
    The tck2connectome command from the mtrix package.

    """
    exe = 'tck2connectome'


class tckgen(BaseMtrixCLI):
    """
    The tckgen command from the mtrix package.

    """
    exe = 'tckgen'


class tckmap(BaseMtrixCLI):
    """
    The tckmap command from the mtrix package.

    """
    exe = 'tckmap'


class tcksift(BaseMtrixCLI):
    """
    The tcksift command from the mtrix package.

    """
    exe = 'tcksift'



def extract_bzero(in_, out):
    return [dwiextract.exe, dwiextract.Flags.bzero, in_, out]
