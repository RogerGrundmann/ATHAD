#!/usr/bin/env python
# ATHAD atmosphere, single epoch, from scratch (restart_from_iter = -1).
from pyathad import Atmosphere

print("=== ATHAD: Hadean atmosphere, from scratch ===", flush=True)
Atmosphere().run()
print("=== ATHAD done ===", flush=True)
