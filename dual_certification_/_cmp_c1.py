from lemma_K import C1_C2_constants, ball_mid
import sys
sys.path.insert(0, "dual_certification_")
from lemma_K import C1_C2_constants, ball_mid
for label, kw in [
  ("conservative", dict(sharp_geom=False, r=6.0, Y0=0.8)),
  ("sharp_default", dict(sharp_geom=True, r=6.0, Y0=0.8)),
  ("sharp_Y0_1.5", dict(sharp_geom=True, r=6.0, Y0=1.5)),
  ("sharp_U2", dict(sharp_geom=True, r=6.0, Y0=1.5, U_norm=2.0)),
]:
  d = C1_C2_constants(field="i", **kw)
  print(label, "C1", ball_mid(d["C1"]), "eta0", ball_mid(d["eta0"]), "C2", ball_mid(d["C2"]), "Ctr", ball_mid(d["C_trace"]))
