"""Out-of-box patpy adapter run on Stephenson: EBMF() with pure defaults.
New defaults (point_laplace, K=10, var_type=1, backfit, seed=67) should
reproduce the tuned result (0.648 disease kNN) exactly.
Usage: python run_stephenson_outofbox_adapter.py <arm_dir> <out_csv>
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
warnings.filterwarnings("ignore")


def _fit_one(item):
    ct, df = item
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
    os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from patpy_cellprograms import EBMF
    m = EBMF(sample_key="sample", cell_group_key="ct")  # pure defaults
    m._fit_in_r({ct: df})
    m._fitted = True
    return ct, m.build_representation()


if __name__ == "__main__":
    arm_dir, out_csv = sys.argv[1], sys.argv[2]
    pb = {f[len("lognorm_hvg__"):-4]: pd.read_csv(os.path.join(arm_dir, f), index_col=0)
          for f in sorted(os.listdir(arm_dir))
          if f.startswith("lognorm_hvg__") and f.endswith(".csv")}
    print(f"cell types: {len(pb)}", flush=True)
    blocks, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=8) as ex:
        for ct, block in ex.map(_fit_one, [(ct, pb[ct]) for ct in sorted(pb)]):
            blocks.append(block)
    rep = pd.concat(blocks, axis=1).fillna(0.0)
    rep.columns = rep.columns.astype(str)
    rep.index = rep.index.astype(str)
    rep.to_csv(out_csv)
    print(f"fit done ({time.time()-t0:.0f}s), rep={rep.shape}", flush=True)
