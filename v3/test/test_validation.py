import unittest
from src.validation import evaluate_model
import numpy as np

class TestValidation(unittest.TestCase):
    def test_evaluate_model(self):
        y_true = np.array([0, 1, 0, 1])
        y_pred = np.array([0, 1, 1, 0])
        results = evaluate_model(y_true, y_pred)
        self.assertTrue('accuracy' in results)
        self.assertTrue('precision' in results)
        self.assertTrue('recall' in results)
        self.assertTrue('f1' in results)


if __name__ == '__main__':
    unittest.main()
