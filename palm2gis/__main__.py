"""Allow `python -m palm2gis` as an alternative to run_palm2gis.py."""
import sys
from .pipeline import main

if __name__ == "__main__":
    sys.exit(main())
