#!/usr/bin/env python3
"""palm2gis entry point.

Usage:
    python run_palm2gis.py -c config/<project>.yaml [-n] [-v|-q]
                           [--stages voxel surfaces wind] [--log-datetime]
"""
import sys
from palm2gis.pipeline import main

if __name__ == "__main__":
    sys.exit(main())
