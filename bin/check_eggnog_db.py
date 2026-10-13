#!/usr/bin/env python
"""Exit non-zero unless the eggNOG DB in <db_dir> is the expected version."""
import sqlite3
import sys

db_dir, expected = sys.argv[1], sys.argv[2]
found = sqlite3.connect(f"{db_dir}/eggnog.db").execute("SELECT * FROM version").fetchone()[0]

if found != expected:
    sys.exit(f"eggNOG DB in {db_dir} is version {found}, pipeline expects {expected}")
print(f"eggNOG DB version {found} OK")
