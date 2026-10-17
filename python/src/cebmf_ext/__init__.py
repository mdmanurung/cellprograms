"""cebmf_ext: scalable cEBMF adapter for cellprograms (workstream G).

Layering rule (binding): this package imports cebmf_torch; nothing imports
this package from cebmf_torch, and cellprograms (R) consumes only the result
contract — never a direct import.
"""
from cebmf_ext.contract import load_result, save_result, validate_contract
from cebmf_ext.fit import FitConfig, FitResult, effective_k, fit_cebmf
from cebmf_ext.preproc import ShiftedLog

__version__ = "0.1.0"
__all__ = [
    "FitConfig",
    "FitResult",
    "ShiftedLog",
    "effective_k",
    "fit_cebmf",
    "load_result",
    "save_result",
    "validate_contract",
]
