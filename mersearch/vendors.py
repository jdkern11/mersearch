from typing import Dict

import selenium

def search_vwr(mol: str) -> Dict[str, str]:
    """Search VWR International website for specific molecule

    Args:
        mol: mol file format string to search for
    Returns:
        Dictionary with data parsed from the VWR search. If no data found, dictionary
        will be empty.
    """

