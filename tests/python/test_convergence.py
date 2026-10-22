import importlib.util
from pathlib import Path
import unittest

path = Path(__file__).resolve().parents[2] / "inst/python/src/cebmf_ext/convergence.py"
spec = importlib.util.spec_from_file_location("convergence", path)
convergence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(convergence)


class ConvergenceTests(unittest.TestCase):
    def test_monotone_plateau_and_flat_loss(self):
        self.assertTrue(convergence.loss_plateau([100, 99, 98.99999, 98.99998], 1e-5, 2))
        self.assertTrue(convergence.loss_plateau([0, 0, 0], 1e-5, 2))

    def test_worsening_never_converges(self):
        self.assertFalse(convergence.loss_plateau([10, 11, 12], 1e-5, 2))
        self.assertFalse(convergence.loss_plateau([10, 10.000001, 10.000002], 1e-5, 2))
        self.assertFalse(convergence.loss_plateau([10, 9.999999, 10], 1e-5, 2))

    def test_nonfinite_short_and_strong_improvement(self):
        for t in ([10], [10, 9, 8], [float("nan"), 1, 0], [10, float("inf"), 1]):
            self.assertFalse(convergence.loss_plateau(t, 1e-5, 2))


if __name__ == "__main__":
    unittest.main()
