"""
traceCB: Trace Cell-type specific eQTL effects across populations using Method of Moments
"""

from .gmm import GMM, GMMtissue
from .ldsc import Run_single_LDSC, Run_cross_LDSC

# Version info
__version__ = "0.1.0"
__author__ = "Luca Jiang"
