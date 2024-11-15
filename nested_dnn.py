import pandas as pd
import sys
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.svm import SVR
from sklearn.metrics import r2_score
from joblib import dump
from sklearn.model_selection import ParameterGrid
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import ExtraTreesRegressor
import time

def nested_cv(data_file, samples):
    # Read the file
    dataset = pd.read_csv(data_file, sep=",", header=None, dtype='str')
    print('I read the file')
    new_dataset1 = dataset.dropna()
    print('I cleaned NAs')
    new_dataset = new_dataset1.sample(frac=1)

    # Prepare hold-out dataset for chromosome 22
    chr22_data_matrix = new_dataset[new_dataset.iloc[:, 0].str.startswith('chr22')]
    new_data_matrix = new_dataset[~new_dataset.iloc[:, 0].str.startswith('chr22')]

    sample_names_chr22 = chr22_data_matrix.iloc[:, 0:1]
    x_chr22 = chr22_data_matrix.iloc[:, 1:-1].values.astype(float).astype(int)
    y_chr22 = chr22_data_matrix.iloc[:, -1]
    unlisted_series_chr22 = y_chr22.apply(lambda x: x[1:-1])
    unlisted_series_chr22 = unlisted_series_chr22.apply(lambda x: '0' if x == 'nan' else x)
    new_y_chr22 = unlisted_series_chr22.values

    # Normalize chr22 data
    Scaler = StandardScaler()
    x_chr22_scaled = Scaler.fit_transform(x_chr22)
    y_chr22_scaled = Scaler.fit_transform(new_y_chr22.reshape(len(new_y_chr22), 1))
    y_chr22_reshaped = np.ravel(y_chr22_scaled)

    # Prepare the training dataset
    sample_names = new_data_matrix.iloc[:, 0:1]
    target = new_data_matrix.iloc[:, -1]
    x_input = new_data_matrix.iloc[:, 1:-1]
    matrix = pd.concat([x_input, target], axis=1).reset_index(drop=True)
    matrix.columns = range(matrix.columns.size)
    X = matrix.iloc[:, :-1].values.astype(float).astype(int)
    y = matrix.iloc[:, -1]
    unlisted_series = y.apply(lambda x: x[1:-1])
    unlisted_series = unlisted_series.apply(lambda x: '0' if x == 'nan' else x)
    new_y = unlisted_series.values

    print('This is the Nested Cross Validation')
    results = []

    # Define the hyperparameters to search for MLPRegressor
    mlp_parameters = {
        'hidden_layer_sizes': [ (360, 360, 360 ), (360, 360,360, 360), (360, 360,360, 360, 360)],
        'activation': ['relu','logistic'],
        'solver': ['adam'],
        'alpha': [0.0001],
        'learning_rate': ['adaptive'],
        'max_iter': [600]
    }

    # Generate all combinations of hyperparameter values for MLPRegressor
    mlp_configurations = list(ParameterGrid(mlp_parameters))

    print('We are searching these parameters: ', mlp_configurations)

    # Perform Nested Cross-Validation (Outer loop with 4 splits)
    final_best_score = -1
    final_best_parameters = None
    outer_loop = KFold(n_splits=4, random_state=None, shuffle=False)
    start_time = time.time()

    for train_index, test_index in outer_loop.split(X):
        print('outer_loop')
        best_score = -1
        best_params = None
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = new_y[train_index], new_y[test_index]

        # Inner loop with 3 splits for hyperparameter search
        inner_cv = KFold(n_splits=3, shuffle=False)
        for configuration in mlp_configurations:
            r2_scores = []
            for train_index_inner, val_index in inner_cv.split(X_train):
                print('inner_loop')
                # Perform train/validation split and normalization
                Scaler = StandardScaler()
                X_train_inner = Scaler.fit_transform(X_train[train_index_inner])
                X_val_inner = Scaler.transform(X_train[val_index])
                y_train_inner = Scaler.fit_transform(y_train[train_index_inner].reshape(-1, 1)).ravel()
                y_val_inner = Scaler.transform(y_train[val_index].reshape(-1, 1)).ravel()

                # Fit the model
                model = MLPRegressor(**configuration)
                model.fit(X_train_inner, np.ravel(y_train_inner))

                # Evaluate the model
                r2 = r2_score(y_val_inner, model.predict(X_val_inner))
                r2_scores.append(r2)
            configuration_score = np.mean(r2_scores)
            if configuration_score > best_score:
                best_score = configuration_score
                best_params = configuration

        print("Best inner accuracy:", best_score, 'Best parameters:', best_params)

        # Outer loop normalization
        Scaler = StandardScaler()
        X_train_normalized = Scaler.fit_transform(X_train)
        X_test_normalized = Scaler.transform(X_test)
        y_train_normalized = Scaler.fit_transform(y_train.reshape(-1, 1)).ravel()
        y_test_normalized = Scaler.transform(y_test.reshape(-1, 1)).ravel()

        # Train the model with the best parameters on the outer training set
        best_model = MLPRegressor(**best_params)
        best_model.fit(X_train_normalized, np.ravel(y_train_normalized))
        test_score = r2_score(y_test_normalized, best_model.predict(X_test_normalized))

        results.append((test_score, best_params))

        if test_score > final_best_score:
            final_best_score = test_score
            final_best_parameters = best_params

        print("Outer fold test accuracy:", final_best_score, 'Final best parameters:', final_best_parameters)

    # Train final model on the entire dataset
    final_best_model = MLPRegressor(**final_best_parameters)
    Scaler = StandardScaler()
    X_normalized = Scaler.fit_transform(X)
    y_normalized = Scaler.fit_transform(new_y.reshape(-1, 1)).ravel()
    final_best_model.fit(X_normalized, np.ravel(y_normalized))

    # Test on chromosome 22 data
    y_chr22_predicted = final_best_model.predict(x_chr22_scaled)
    r2 = r2_score(y_chr22_reshaped, y_chr22_predicted)
    print('Chromosome 22 R² score:', r2)
    print("Final test accuracy:", final_best_score, 'Final best parameters:', final_best_parameters)
    
    # Save the best model
    dump(final_best_model, 'best_model_nested_cv_dnn_h3k4me3' + str(samples) + '.joblib')
    end_time = time.time()
    print('Total time:', end_time - start_time)
    return "Done!"
cross_val_dataset = sys.argv[1]
samples =  sys.argv[2]
#hold_out_dataset =  sys.argv[3]
print(nested_cv(cross_val_dataset, samples))