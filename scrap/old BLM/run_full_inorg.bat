blm_245d.exe /S full_inorg_TOX.blmscr
del full_inorg_TOX.det.xls
del full_inorg_TOX.sim.xls
ren full_inorg.det.xls full_inorg_TOX.det.xls
ren full_inorg.sim.xls full_inorg_TOX.sim.xls

blm_245d.exe /S full_inorg_SPEC.blmscr
del full_inorg_SPEC.det.xls
del full_inorg_SPEC.sim.xls
ren full_inorg.det.xls full_inorg_SPEC.det.xls
ren full_inorg.sim.xls full_inorg_SPEC.sim.xls

pause