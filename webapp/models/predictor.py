import pandas as pd
import sklearn
import pickle
import numpy as np


def get_predictions(entry,):
    try:
        with open("/var/www/tmsnp/tmsnp/models/estimator.pkl", "rb") as input_file:
            dict_estimator = pickle.load(input_file)
    except FileNotFoundError:
        raise FileNotFoundError

    estimator = dict_estimator["estimator"]
    scaler = dict_estimator["scaler"]
    entry = scaler.transform(entry)
    prediction = estimator.predict(entry, 0.2)

    return prediction


# test = pd.read_csv('../../test_ref.txt', sep='\t')
# entry = np.asarray([test.iloc[0][:-1]])
# print(entry)
# pred = get_predictions(entry)
# print(pred)
