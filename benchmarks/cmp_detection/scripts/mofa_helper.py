"""mofa_helper.py -- minimal MOFA+ (mofapy2) wrapper for the CMP v2 joint-factorization arm.

fit_mofa(views, K, seed, maxiter) -> dict(W=[features x K per view], Z=samples x K, r2=K x V)
views: list of 2D float arrays (samples x features), identical sample order, NaN rows for missing samples.
"""
import os, sys, io, contextlib
import numpy as np


def fit_mofa(views, K=15, seed=1, maxiter=1000):
    from mofapy2.run.entry_point import entry_point
    views = [np.asarray(v, dtype=float) for v in views]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ent = entry_point()
        ent.set_data_options(scale_views=True, scale_groups=False, center_groups=True)
        ent.set_data_matrix([[v] for v in views], likelihoods=["gaussian"] * len(views))
        ent.set_model_options(factors=int(K), spikeslab_weights=True, ard_weights=True, ard_factors=False)
        ent.set_train_options(iter=int(maxiter), convergence_mode="fast", startELBO=1, freqELBO=5,
                              dropR2=None, gpu_mode=False, seed=int(seed), verbose=False, quiet=True)
        ent.build()
        ent.run()
        W = [np.asarray(w) for w in ent.model.nodes["W"].getExpectation()]
        Z = np.asarray(ent.model.nodes["Z"].getExpectation())
        Kf = Z.shape[1]
        # per-view variance explained by each factor (centered data, missing ignored)
        r2 = np.zeros((Kf, len(views)))
        Ys = ent.data
        for m in range(len(views)):
            Y = np.asarray(Ys[m]); ok = ~np.isnan(Y).any(axis=1)
            Y = Y[ok]; Zm = Z[ok]
            ss = np.nansum(Y ** 2)
            for k in range(Kf):
                pred = np.outer(Zm[:, k], W[m][:, k])
                r2[k, m] = 1.0 - np.nansum((Y - pred) ** 2) / ss
    return {"W": W, "Z": Z, "r2": r2}
