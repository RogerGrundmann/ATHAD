#!/usr/bin/env python
# Ma=0 atmosphere, from scratch (restart_from_iter=-1), to iter 400. Writes Transfer_Atm_<iter>.
from pyatom import Atmosphere
print("=== ATM Ma=0 from scratch -> iter 400 ===", flush=True)
Atmosphere().run()
print("=== ATM done ===", flush=True)
