### Perform cross validation for other kind od regressors such as extra trees or random forest
from sklearn.feature_selection import SelectKBest
from sklearn.datasets import make_regression
from sklearn.feature_selection import r_regression
import pandas as pd
import sys
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import  r2_score, mean_squared_error, mean_absolute_error
from sklearn.model_selection import cross_validate
from sklearn.model_selection import GridSearchCV
from joblib import dump, load
import time
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.model_selection import KFold
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import ParameterGrid
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import RandomForestRegressor


def cv_hold_out(data_file, model_filename):
    # read the file
    dataset = pd.read_csv(data_file, sep=",", header=None, dtype='str')
    print('I read the file')
    new_dataset1=dataset.dropna()
    print('I cleaned NAs')
    new_dataset = new_dataset1.sample(frac=1)
    
    # Prepare hold out dataset , for testing of chrosome 22 values
    chr7_data_matrix = new_dataset[new_dataset.iloc[:, 0].str.startswith('chr7')]
    new_data_matrix = new_dataset[~new_dataset.iloc[:, 0].str.startswith('chr7')]
    sample_names_chr7 = chr7_data_matrix.iloc[:, 0]
    x_chr7 = chr7_data_matrix.iloc[:,1:-1].values.astype(float).astype(int)
    y_chr7 =  chr7_data_matrix.iloc[:,-1]
    unlisted_series_chr7 = y_chr7.apply(lambda x: x[1:-1])
    unlisted_series_chr7 = unlisted_series_chr7.apply(lambda x: '0' if x == 'nan' else x)
    print(unlisted_series_chr7)
    new_y_chr7 = unlisted_series_chr7.values

    Scaler = StandardScaler()
    x_chr7_scaled=Scaler.fit_transform(x_chr7[:,:-3])
    y_chr7_scaled = Scaler.fit_transform(new_y_chr7.reshape(len(new_y_chr7), 1))
    y_chr7_reshaped = np.ravel(y_chr7_scaled)
    whole_x_chr7 = (np.column_stack((x_chr7_scaled, x_chr7[:,-3:])))
    
    # Prepare training matrix
    sample_names = new_data_matrix.iloc[:, 0]
    sample_names = sample_names.str.split('_', expand=True)
    print(sample_names)
    target=new_dataset.iloc[:,-1]
    x_input = new_dataset.iloc[:,1:-1]
    matrix=pd.concat([x_input,sample_names.iloc[:,3:],target],axis=1).reset_index(drop=True)
    print(matrix)
    matrix.columns=range(matrix.columns.size)
    x = matrix.iloc[:,:-1].values.astype(float).astype(int)
    y = matrix.iloc[:,-1]
    unlisted_series = y.apply(lambda x: x[1:-1])
    unlisted_series = unlisted_series.apply(lambda x: '0' if x == 'nan' else x)
    print(unlisted_series)
    new_y = unlisted_series.values
    
    X_train, X_test, y_train, y_test = train_test_split(x, new_y, test_size=0.25, random_state=42)
    model_trees = ExtraTreesRegressor(n_estimators=600, criterion='squared_error', min_samples_split = 2,min_samples_leaf = 1, max_depth = None)
    dnn = MLPRegressor(random_state=42, max_iter=400, hidden_layer_sizes=(362, 362, 362, 362))
    Scaler = StandardScaler()
    print(X_train[:,:-3])
    x_scaled, x_test_scaled = Scaler.fit_transform(X_train[:,:-3]), Scaler.fit_transform(X_test[:,:-3])
    y_scaled, y_test_scaled = Scaler.fit_transform(y_train.reshape(len(y_train), 1)), Scaler.fit_transform(y_test.reshape(len(y_test), 1))
    y_reshaped = np.ravel(y_scaled)
    y_reshaped_test = np.ravel(y_test_scaled)
    print(np.column_stack((x_scaled, X_train[:,-3:])))
    model_trees.fit(np.column_stack((x_scaled, X_train[:,-3:])), y_reshaped)
    trees_score = model_trees.score(np.column_stack((x_test_scaled, X_test[:,-3:])), y_reshaped_test)
    dnn.fit(np.column_stack((x_scaled, X_train[:,-3:])), y_reshaped)
    dnn_score = dnn.score(np.column_stack((x_test_scaled, X_test[:,-3:])), y_reshaped_test)
    print ("Score is 500 trees:", trees_score , 'dnn 4*362, 400:',  dnn_score)
    dump(dnn, 'dnn_model_h2bub2.joblib')
    dump(model_trees, model_filename)


    trees_scoreChr7 = model_trees.score(whole_x_chr7, y_chr7_reshaped)
    dnn_scoreChr7 = dnn.score(whole_x_chr7, y_chr7_reshaped)
    print(trees_scoreChr7, dnn_scoreChr7)


def sum_results(results_sorted):
    steps = [ 1000, 100]
    for step in steps:
        Y_original = []
        Y_predicted = []
        print(step)
        print(len(Y_original))
        for i in range(int(len(results_sorted)/step)):
            #print(step*i,(step*i)+step)
            Y_original.append(sum(results_sorted.iloc[(step*i):(step*i)+step,]['Y_original_scaled']))
            Y_predicted.append(sum(results_sorted.iloc[(step*i):(step*i)+step,]['Y_predicted']))
        print(r2_score(Y_original, Y_predicted))
    r2_last_1kb = r2_score(Y_original, Y_predicted)
    return(r2_last_1kb)

def model_another_dataset(data_matrix_condition, model_filename, output_filename):
    dataset = pd.read_csv(data_matrix_condition, sep=",", header=None, dtype='str')
    print('I read the file')
    new_dataset1=dataset.dropna()
    print('I cleaned NAs')
    new_dataset = new_dataset1.sample(frac=1)
    # Predict for other dataset that corresponds to another condition (heatshock for now)
    sample_names_condition = new_dataset.iloc[:, 0]
    sample_names = sample_names_condition.str.split('_', expand=True)
    print(sample_names)
    target=new_dataset.iloc[:,-1]
    x_input = new_dataset.iloc[:,1:-1]
    matrix=pd.concat([x_input,sample_names.iloc[:,3:],target],axis=1).reset_index(drop=True)
    print(matrix)
    matrix.columns=range(matrix.columns.size)
    x_co = matrix.iloc[:,:-1].values.astype(float).astype(int)
    y_co = matrix.iloc[:,-1]
    unlisted_series = y_co.apply(lambda x: x[1:-1])
    unlisted_series = unlisted_series.apply(lambda x: '0' if x == 'nan' else x)
    print(unlisted_series)
    y_co = unlisted_series.values

    #sample_names_condition = new_dataset.iloc[:, 0]
    #target=new_dataset.iloc[:,-1]
    #x_input = new_dataset.iloc[:,1:-1]
    #matrix=pd.concat([x_input,target],axis=1).reset_index(drop=True)
    #print(matrix)
    #matrix.columns=range(matrix.columns.size)
    #x_co = matrix.iloc[:,:-1].values.astype(float).astype(int)
    #y_co = matrix.iloc[:,-1]
    #unlisted_series =(y_co.apply(lambda x: x[1:-1]))
    #y_without_nan = unlisted_series.apply(lambda x: '0' if x == 'nan' else x)
    #print(y_without_nan)
    #y_co = y_without_nan.values

    Scaler_co = StandardScaler()
    x_scaled_co=Scaler_co.fit_transform(x_co)
    y_scaled_co=Scaler_co.fit_transform(y_co.reshape(len(y_co), 1))
    regressor = load(model_filename)
    y_predicted_co = regressor.predict(x_scaled_co)
    r2 = r2_score(y_scaled_co, y_predicted_co)
    mae = mean_absolute_error(y_scaled_co, y_predicted_co)
    mse = mean_squared_error(y_scaled_co, y_predicted_co)
    print('condition scores of the predict: ', r2, mae, mse)
    data_co = {
    'samples': (sample_names_condition),
    'Y_original_scaled': list(np.ravel(y_scaled_co)),
    'Y_predicted':  list(y_predicted_co)}
    df_co = pd.DataFrame(data_co)
    df_co.to_csv(output_filename, sep='\t', index=False)

    #results_sorted = df_co.sort_values(by = ['start'])
    #df_filtered = results_sorted[results_sorted['Y_original_scaled'] != 0.0]
    #rint(sum_results(df_filtered))
    #print(sum_results(results_sorted))

wildtype_dataset = sys.argv[1]
condition_dataset = sys.argv[2]
model_filename = sys.argv[3]
output_filename = sys.argv[4]

#print(cv_hold_out(wildtype_dataset, model_filename))
print(model_another_dataset(condition_dataset, model_filename, output_filename))
