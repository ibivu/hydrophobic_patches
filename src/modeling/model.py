from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
import numpy as np
import xgboost
import pickle
import yaml
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
import pandas as pd

def get_best_model(X, y, methods, model_output=None):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    highest_r = 0
    best_model = None

    for method in methods:
        clf = GridSearchCV(methods[method][0], methods[method][1], cv=5)
        clf.fit(X_train, y_train)
        rsquared = clf.score(X_test, y_test)

        if rsquared > highest_r:
            highest_r = rsquared
            best_model = clf

    best_model = best_model.fit(X,y)
    if model_output:
        pickle.dump(best_model, model_output)


    return best_model



config = yaml.safe_load(open("../config.yml"))

svr = SVR(kernel='rbf')
svr_params = {"C": [1e0, 1e1, 1e2, 1e3], "gamma": np.logspace(-2, 2, 5)}

xgbbooster = xgboost.XGBRegressor(objective='reg:squarederror')
xgbbooster_params = {
        'min_child_weight': [1, 5, 10],
        'gamma': [0, 0.5, 1.5, 5],
        # 'subsample': [0.6, 0.8, 1.0],
        # 'colsample_bytree': [0.6, 0.8, 1.0],
        # 'max_depth': [2, 3, 4]
        }

randfor = RandomForestRegressor()
randfor_params = {'max_depth':[2,3,4]}

destree = DecisionTreeRegressor()
destree_params = {'max_depth':[2,3,4]}

linreg = LinearRegression()
linreg_params = {}

methods = {
'xgboost':[xgbbooster, xgbbooster_params],
'svr':[svr, svr_params],
'destree':[destree, destree_params],
'randfor':[randfor,randfor_params],
}

df_train = pd.read_csv(config['path']['processed_data']+'ready_to_use_data_train.csv')
df_test = pd.read_csv(config['path']['processed_data']+'ready_to_use_data_test.csv')

gfm_columns = [ 'length', 'entropy',
       'hydr_count', 'polar_count', 'buried',
       'gravy', 'molecular_weight', 'aromaticity', 'instability_index',
       'isoelectric_point', 'A', 'C', 'D', 'E', 'F',
       'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
       'Y']
length = ['length']
nsp2_columns = ['rhsa_netsurfp2','thsa_netsurfp2','tasa_netsurfp2']
nsp2_gfm_columns =  nsp2_columns + gfm_columns

X_train_gfm = df_train[gfm_columns]
X_train_nsp2 = df_train[nsp2_columns]
X_train_nsp2_gfm = df_train[nsp2_gfm_columns]
X_train_length = df_train[length]
y_train_thsa, y_train_rhsa, y_train_lhpsa = df_train['thsa'], df_train['rhsa'], df_train['size']

X_test_gfm = df_test[gfm_columns]
X_test_nsp2 = df_test[nsp2_columns]
X_test_nsp2_gfm = df_test[nsp2_gfm_columns]
X_test_length = df_test[length]
y_test_thsa, y_test_rhsa, y_test_lhpsa = df_test['thsa'], df_test['rhsa'], df_test['size']
y_test_id = df_test['id']

thsa_gfm_model = get_best_model(X_train_gfm, y_train_thsa, methods, open(config['path']['model']+'thsa_gfm.model', 'wb'))
thsa_nsp2_gfm_model = get_best_model(X_train_nsp2_gfm, y_train_thsa, methods, open(config['path']['model']+'thsa_nsp2_gfm.model', 'wb'))
thsa_length_model = get_best_model(X_train_length, y_train_thsa, {'linreg':[linreg, linreg_params]}, open(config['path']['model']+'thsa_length.model', 'wb'))

rhsa_gfm_model = get_best_model(X_train_gfm, y_train_rhsa, methods, open(config['path']['model']+'rhsa_gfm.model', 'wb'))
rhsa_nsp2_gfm_model = get_best_model(X_train_nsp2_gfm, y_train_rhsa, methods, open(config['path']['model']+'rhsa_nsp2_gfm.model', 'wb'))
rhsa_length_model = get_best_model(X_train_length, y_train_rhsa, {'linreg':[linreg, linreg_params]}, open(config['path']['model']+'rhsa_length.model', 'wb'))

lhpsa_gfm_model = get_best_model(X_train_gfm, y_train_lhpsa, methods, open(config['path']['model']+'lhpsa_gfm.model', 'wb'))
lhpsa_nsp2_model = get_best_model(X_train_nsp2, y_train_lhpsa, methods, open(config['path']['model']+'lhpsa_nsp2.model', 'wb'))
lhpsa_nsp2_gfm_model = get_best_model(X_train_nsp2_gfm, y_train_lhpsa, methods, open(config['path']['model']+'lhpsa_nsp2_gfm.model', 'wb'))
lhpsa_length_model = get_best_model(X_train_length, y_train_lhpsa, {'linreg':[linreg, linreg_params]}, open(config['path']['model']+'lhpsa_length.model', 'wb'))

pd.DataFrame({'id':y_test_id, 'prediction':thsa_gfm_model.predict(X_test_gfm)}).to_csv(config['path']['predictions']+'thsa_gfm_prediction.csv', index=False)
pd.DataFrame({'id':y_test_id, 'prediction':thsa_nsp2_gfm_model.predict(X_test_nsp2_gfm)}).to_csv(config['path']['predictions']+'thsa_nsp2_gfm_prediction.csv', index=False)
pd.DataFrame({'id':y_test_id, 'prediction':thsa_length_model.predict(X_test_length)}).to_csv(config['path']['predictions']+'thsa_length_prediction.csv', index=False)

pd.DataFrame({'id':y_test_id, 'prediction':rhsa_gfm_model.predict(X_test_gfm)}).to_csv(config['path']['predictions']+'rhsa_gfm_prediction.csv', index=False)
pd.DataFrame({'id':y_test_id, 'prediction':rhsa_nsp2_gfm_model.predict(X_test_nsp2_gfm)}).to_csv(config['path']['predictions']+'rhsa_nsp2_gfm_prediction.csv', index=False)
pd.DataFrame({'id':y_test_id, 'prediction':rhsa_length_model.predict(X_test_length)}).to_csv(config['path']['predictions']+'rhsa_length_prediction.csv', index=False)

pd.DataFrame({'id':y_test_id, 'prediction':lhpsa_gfm_model.predict(X_test_gfm)}).to_csv(config['path']['predictions']+'lhpsa_gfm_prediction.csv', index=False)
pd.DataFrame({'id':y_test_id, 'prediction':lhpsa_nsp2_model.predict(X_test_nsp2)}).to_csv(config['path']['predictions']+'lhpsa_nsp2_prediction.csv', index=False)
pd.DataFrame({'id':y_test_id, 'prediction':lhpsa_nsp2_gfm_model.predict(X_test_nsp2_gfm)}).to_csv(config['path']['predictions']+'lhpsa_nsp2_gfm_prediction.csv', index=False)
pd.DataFrame({'id':y_test_id, 'prediction':lhpsa_length_model.predict(X_test_length)}).to_csv(config['path']['predictions']+'lhpsa_length_prediction.csv', index=False)
