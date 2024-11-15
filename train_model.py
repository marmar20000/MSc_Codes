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
from thundersvm import SVR as fasterSVR
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.model_selection import cross_val_score
from sklearn.metrics import  r2_score
import pandas as pd
import os
import numpy as np
import math
import sys
import random


def cv_hold_out(data_file, model_filename):
    # read the file
    dataset = pd.read_csv(data_file, sep=",", header=None, dtype='str')
    print('I read the file')
    new_dataset1=dataset.dropna()
    print('I cleaned NAs')
    new_dataset = new_dataset1.sample(frac=1)

    # Prepare hold out dataset , for testing of chrosome 22 values
    chr22_data_matrix = new_dataset[new_dataset.iloc[:, 0].str.startswith('chr22')]
    new_data_matrix = new_dataset[~new_dataset.iloc[:, 0].str.startswith('chr22')]
    sample_names_chr22 = chr22_data_matrix.iloc[:, 0:1]
    x_chr22 = chr22_data_matrix.iloc[:,1:-1].values.astype(float).astype(int)
    y_chr22 =  chr22_data_matrix.iloc[:,-1]
    unlisted_series_chr22 = y_chr22.apply(lambda x: x[1:-1])
    unlisted_series_chr22 = unlisted_series_chr22.apply(lambda x: '0' if x == 'nan' else x)
    print(unlisted_series_chr22)
    new_y_chr22 = unlisted_series_chr22.values
    #new_y_chr22 = y_chr22.values
    Scaler = StandardScaler()
    x_chr22_scaled=Scaler.fit_transform(x_chr22)
    y_chr22_scaled = Scaler.fit_transform(new_y_chr22.reshape(len(new_y_chr22), 1))
    y_chr22_reshaped = np.ravel(y_chr22_scaled)

    # Prepare training matrix
    sample_names = new_data_matrix.iloc[:, 0:1]
    target=new_data_matrix.iloc[:,-1]
    x_input = new_data_matrix.iloc[:,1:-1]
    matrix=pd.concat([x_input,target],axis=1).reset_index(drop=True)
    print(matrix)
    matrix.columns=range(matrix.columns.size)
    x = matrix.iloc[:,:-1].values.astype(float).astype(int)
    y = matrix.iloc[:,-1]
    unlisted_series = y.apply(lambda x: x[1:-1])
    unlisted_series = unlisted_series.apply(lambda x: '0' if x == 'nan' else x)
    print(unlisted_series)
    new_y = unlisted_series.values
    #new_y = y.values
    #model_trees = ExtraTreesRegressor(n_estimators=700, criterion='squared_error', n_jobs = -2, min_samples_split = 2, min_impurity_decrease = 0)
    regr = MLPRegressor(random_state=42, alpha=0.0001, learning_rate = 'adaptive',   max_iter=600, hidden_layer_sizes=(360, 360, 360), activation='relu', solver='adam')
    #regr = MLPRegressor(random_state=42, alpha=0.0001, learning_rate = 'adaptive',   max_iter=600, hidden_layer_sizes=(360, 360, 360, 360), activation='relu', solver='adam')
    Scaler = StandardScaler()
    x_scaled = Scaler.fit_transform(x)
    y_scaled = Scaler.fit_transform(new_y.reshape(len(new_y), 1))
    y_reshaped = np.ravel(y_scaled)
    regr.fit(x_scaled, y_reshaped)
    score_chr22_dnn = regr.score(x_chr22_scaled, y_chr22_reshaped)
    #model_trees.fit(x_scaled, y_reshaped)
    #score_chr22_trees = model_trees.score(x_chr22_scaled, y_chr22_reshaped)
    print('the score of the hold out chr 22 with model trained with all dataset is: ',  score_chr22_dnn)
    dump(regr, model_filename)
    print('winner is DNN')


def sum_results(results_sorted):
    steps = [ 1000, 100]
    for step in steps:
        Y_original = []
        Y_predicted = []

        for i in range(int(len(results_sorted)/step)):
            #print(step*i,(step*i)+step)
            Y_original.append(sum(results_sorted.iloc[(step*i):(step*i)+int(step),]['Y_original_scaled']))
            Y_predicted.append(sum(results_sorted.iloc[(step*i):(step*i)+int(step),]['Y_predicted']))
        print(r2_score(Y_original, Y_predicted))
    r2 = r2_score(Y_original, Y_predicted)
    return(r2)


def model_another_dataset(data_matrix_condition, model_filename, output_filename):
    dataset = pd.read_csv(data_matrix_condition, sep=",", header=None, dtype='str')
    print('I read the file')
    new_dataset1=dataset.dropna()
    print('I cleaned NAs')
    new_dataset = new_dataset1.sample(frac=1)
    # Predict for other dataset that corresponds to another condition (heatshock for now)
    sample_names_condition = new_dataset.iloc[:, 0:3]
    target=new_dataset.iloc[:,-1]
    x_input = new_dataset.iloc[:,3:-1]
    matrix=pd.concat([x_input,target],axis=1).reset_index(drop=True)
    print(matrix)
    matrix.columns=range(matrix.columns.size)
    x_co = matrix.iloc[:,:-1].values.astype(float).astype(int)
    y_co = matrix.iloc[:,-1].values
    #unlisted_series =(y_co.apply(lambda x: x[1:-1]))
    #y_without_nan = unlisted_series.apply(lambda x: '0' if x == 'nan' else x)
    #print(y_without_nan)
    #y_co = y_without_nan.values

    Scaler_co = StandardScaler()
    x_scaled_co=Scaler_co.fit_transform(x_co)
    y_scaled_co=Scaler_co.fit_transform(y_co.reshape(len(y_co), 1))
    #print(Scaler_co.mean_, Scaler_co.scale_)
    regressor = load(model_filename)
    y_predicted_co = regressor.predict(x_scaled_co)
    r2 = r2_score(y_scaled_co, y_predicted_co)
    mae = mean_absolute_error(y_scaled_co, y_predicted_co)
    mse = mean_squared_error(y_scaled_co, y_predicted_co)
    print('condition scores of the predict: ', r2, mae, mse)
    data_co = {
    #'Samples' : np.ravel(sample_names_condition   ),
    'chr': (sample_names_condition.iloc[:, 0]),
    'start':(sample_names_condition.iloc[:, 1]),
    'end' : (sample_names_condition.iloc[:, 2]),
    'Y_original_scaled': list(np.ravel(y_scaled_co)),
    'Y_predicted':  list(y_predicted_co)}
    df_co = pd.DataFrame(data_co)
    df_co.to_csv(output_filename+".tsv", sep='\t', index=False)
    dataset = pd.read_csv(output_filename+".tsv", sep = '\t')
    dataset_sorted = dataset.sort_values(by = ['start'])
    df_filtered = dataset_sorted[dataset_sorted['Y_original_scaled'] != 0.0]
    print(sum_results(dataset_sorted))
    print(sum_results(df_filtered))



def only_evaluation(file):
    dataset = pd.read_csv(file, sep = '\t')
    dataset_sorted = dataset.sort_values(by = ['start'])
    df_filtered = dataset_sorted[dataset_sorted['Y_original_scaled'] != 0.0]
    print(sum_results(dataset_sorted))
    print(sum_results(df_filtered))



def bed_to_wig(bed_file, output_file):
    data = {}  # Dictionary to store aggregated scores
    # Read BED file
    with open(bed_file, 'r') as f:
      
        for line in f:
        
            chrom, start, end, experiment, score = line.split('\t')
            start, end, score = int(start), int(end), float(score)
            for pos in range(start, end):
                if chrom not in data:
                    data[chrom] = {}
                if pos not in data[chrom]:
                    data[chrom][pos] = []
                data[chrom][pos].append(score)
    # Write data to WIG file
    with open(output_file, 'w') as f:
        f.write("track type=wiggle_0 name="+output_file+" description=\"BED to WIG conversion\"\n")
        for chrom in data:
            f.write("variableStep chrom={}\n".format(chrom))
            for pos in sorted(data[chrom]):
                score = max(data[chrom][pos])  # Aggregate score (e.g., taking max)
                f.write("{} {}\n".format(pos+1, score))  # Position is 1-based
    

def ready_for_plotting(file_tomanilupate, predicted_out_plot, original_out_plot):
    with open(file_tomanilupate , 'r') as file:
        y_predicted =[]
        y_original_scaled = []
        with open(predicted_out_plot, 'w') as bed_predicted:
            for i in file.readlines()[1:]:
                y_predicted.append((float((i.strip('\n').split('\t'))[2])))
                y_original_scaled.append(float((i.strip('\n').split('\t'))[1]))
                bed_predicted.write(i.strip('\n').split('\t')[0].split("_")[0]+'\t'+i.strip('\n').split('\t')[0].split("_")[1]+'\t'+i.strip('\n').split('\t')[0].split("_")[2]+'\t'+"predicted"+'\t'+str(float(i.strip('\n').split('\t')[2]))+'\n')
    r2 = r2_score(y_original_scaled, y_predicted)
    print(r2)
    with open(file_tomanilupate , 'r') as file:
        y_predicted =[]
        y_original_scaled = []
        with open( original_out_plot, 'w') as bed_original:
            for i in file.readlines()[1:]:
                bed_original.write(i.strip('\n').split('\t')[0].split("_")[0]+'\t'+i.strip('\n').split('\t')[0].split("_")[1]+'\t'+i.strip('\n').split('\t')[0].split("_")[2]+'\t'+"original"+'\t'+str(float(i.strip('\n').split('\t')[1]))+'\n')


wildtype_dataset = sys.argv[1]
condition_dataset = sys.argv[2]
model_filename = sys.argv[3]
output_filename = sys.argv[4]

print(cv_hold_out(wildtype_dataset, model_filename))
print(model_another_dataset(condition_dataset, model_filename, output_filename))
ready_for_plotting(output_filename+".tsv", output_filename+"_predicted.bed", output_filename+"_original.bed")
bed_to_wig(output_filename+"_predicted.bed", output_filename+"_predicted.wig")
bed_to_wig(output_filename+"_original.bed", output_filename+"_original.wig")
