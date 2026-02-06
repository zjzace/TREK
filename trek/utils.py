"""
Utility functions for polyAFinder
(Currently not used in the pipeline - reserved for future extensions)
"""

import logging

logger = logging.getLogger(__name__)


def format_position(pos: int, to_one_based: bool = True) -> int:
    """
    Convert between 0-based and 1-based coordinates
    
    Args:
        pos: Position
        to_one_based: If True, convert to 1-based; if False, convert to 0-based
        
    Returns:
        Converted position
    """
    return pos + 1 if to_one_based else pos - 1

