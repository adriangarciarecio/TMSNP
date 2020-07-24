import pandas as pd
import sklearn
import pickle
import numpy as np


def get_predictions(entry, conf, mode):
    try:
        if mode == "pathogenic":
            with open(
                "/var/www/tmsnp/tmsnp/models/estimator_patho_v2.pkl", "rb"
            ) as input_file:
                dict_estimator = pickle.load(input_file)
        if mode == "nonpathogenic":
            with open(
                "/var/www/tmsnp/tmsnp/models/estimator_nonpatho_v2.pkl", "rb"
            ) as input_file:
                dict_estimator = pickle.load(input_file)
    except FileNotFoundError:
        raise FileNotFoundError

    estimator = dict_estimator["estimator"]
    scaler = dict_estimator["scaler"]
    entry = scaler.transform(entry)
    prediction = estimator.predict(entry, conf)
    if prediction.item((0, 0)) == True:
        if prediction.item((0, 1)) == False:
            pred = "Positive"
        else:
            pred = "Outside the domain of applicability"
    if prediction.item((0, 0)) == False:
        if prediction.item((0, 1)) == True:
            pred = "Negative"
        else:
            pred = "Outside the domain of applicability"
    return pred


# test = pd.read_csv('../../test_ref.txt', sep='\t')
# entry = np.asarray([test.iloc[0][:-1]])
# print(entry)
# pred = get_predictions(entry)
# print(pred)
