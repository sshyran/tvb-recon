class Atlas(object):
    DEFAULT = "default"
    A2009S = "a2009s"
    DKT = "DKT"


# TODO: better set properties to the exact atlases names to avoid if statements and easily add more parcellations
# i.e.:
# class AtlasSuffixDict(object):
#     default = ""
#     a2009s = ".a2009s"
#     DKT = ".DKTatlas"
# with desired use: atlas_suffix = getattr(AtlasSuffix, current_atlas, None)
class AtlasSuffix(object):
    DEFAULT = ""
    A2009S = ".a2009s"
    DKT = ".DKTatlas"