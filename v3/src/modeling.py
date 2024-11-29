from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

def train_model(X_train, y_train):
    """Trains a random forest classifier."""
    model = RandomForestClassifier()
    model.fit(X_train, y_train)
    return model

def predict(model, X_test):
    """Makes predictions using the trained model."""
    return model.predict(X_test)
