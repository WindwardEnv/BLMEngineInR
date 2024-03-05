blm_245c.exe /S full_organic_TOX.blmscr
del full_organic_TOX.det.xls
del full_organic_TOX.sim.xls
ren full_organic.det.xls full_organic_TOX.det.xls
ren full_organic.sim.xls full_organic_TOX.sim.xls

blm_245c.exe /S full_organic_SPEC.blmscr
del full_organic_SPEC.det.xls
del full_organic_SPEC.sim.xls
ren full_organic.det.xls full_organic_SPEC.det.xls
ren full_organic.sim.xls full_organic_SPEC.sim.xls

pause