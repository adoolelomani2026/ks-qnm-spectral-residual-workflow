#!/usr/bin/env python3
"""Run the Leaver-validated QNM catalogue generator from the repo root."""

from __future__ import annotations

import sys
from pathlib import Path
import os

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
os.chdir(ROOT_DIR)

from qnm.catalogue import main


if __name__ == "__main__":
    main()
