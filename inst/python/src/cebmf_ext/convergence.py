"""Loss-plateau oracle, independent of optional GPU/model dependencies."""
import math


def loss_plateau(trace, tolerance, patience):
    """Only finite, non-increasing small loss changes establish convergence."""
    if patience < 1 or tolerance <= 0 or len(trace) < patience + 1:
        return False
    tail = trace[-(patience + 1):]
    if not all(math.isfinite(v) for v in tail):
        return False
    scale = max(abs(trace[-1]), abs(trace[0]), 1e-300)
    gains = [(a - b) / scale for a, b in zip(tail, tail[1:])]
    return all(0 <= g < tolerance for g in gains)
