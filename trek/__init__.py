"""
polyAFinder: A tool for identifying alternative polyA sites from long-read data
"""

__version__ = "1.0.0"
__author__ = "Your Name"

from .gtf_processor import GTFProcessor
from .alignment_processor import AlignmentProcessor
from .apa_finder import TESFinder

__all__ = ['GTFProcessor', 'AlignmentProcessor', 'TESFinder']
