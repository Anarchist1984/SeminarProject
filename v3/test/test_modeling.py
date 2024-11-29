import unittest
from src.modeling import train_model, predict
import numpy as np
from sklearn.exceptions import NotFittedError

class TestModeling(unittest.TestCase):
    def test_train_model(self):
        X_train = np.array([[1, 2], [3, 4]])
        y_train = np.array([0, 1])
        model = train_model(X_train, y_train)
        self.assertTrue(hasattr(model, 'predict'))

    def test_predict(self):
        X_train = np.array([[1, 2], [3, 4]])
        y_train = np.array([0, 1])
        model = train_model(X_train, y_train)
        X_test = np.array([[5, 6]])
        y_pred = predict(model, X_test)
        self.assertEqual(len(y_pred), 1)
        with self.assertRaises(NotFittedError):
            model = train_model(X_train, y_train)
            predict(model, X_test)

if __name__ == '__main__':
    unittest.main()
