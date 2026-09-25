"""
traceCB: Trace Cell-type specific eQTL effects across populations using Method of Moments
"""

from .gmm import GMM, GMMtissue
from .ldsc import Run_Cross_LDSC
from . import utils

# Version info
__version__ = "1.0"
__author__ = "Wenxin Jiang, Jiashun Xiao, and Mingxuan Cai"

__all__ = [
    "GMM",
    "GMMtissue",
    "Run_Cross_LDSC",
    "utils",
]
