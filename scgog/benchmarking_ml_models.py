# Import necessary libraries
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
import pandas as pd
import muon as mu

from data_loader import get_loader

def benchmark_models(mdata, test_size = 0.2, random_state=42):
    """
    Benchmarks several ML models on a dataset with clustering-based labels.

    Parameters:
    - mdata_path (str): Path to the .h5mu file containing omics data.
    """
    # Load the preprocessed mdata
    
    
    # Load the dataset
    loader = get_loader(mdata)

    # Convert DataLoader to numpy arrays for scikit-learn
    features, labels = next(iter(loader))
    features = features.numpy()
    labels = labels.numpy()

    # Split the dataset
    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=test_size, random_state=random_state)

    # Define models to benchmark
    models = {
        "Random Forest": RandomForestClassifier(n_estimators=100),
        "Gradient Boosting": GradientBoostingClassifier(n_estimators=100),
        "SVC": SVC(gamma='auto'),
        "Logistic Regression": LogisticRegression(max_iter=200)
    }

    # Results container
    results = []

    # Train and evaluate each model
    for name, model in models.items():
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        # Calculate metrics
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred, average='weighted')
        recall = recall_score(y_test, y_pred, average='weighted')
        f1 = f1_score(y_test, y_pred, average='weighted')

        results.append({"Model": name, "Accuracy": accuracy, "Precision": precision, "Recall": recall, "F1 Score": f1})

    # Convert results to DataFrame for nicer display
    results_df = pd.DataFrame(results)
    print(results_df)
