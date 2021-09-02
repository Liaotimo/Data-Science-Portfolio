# Iris-Species-by-Measurements
In this project, the famous [Iris flower data set](http://en.wikipedia.org/wiki/Iris_flower_data_set) will be used to train a model for prediction of iris species. This data set is composed of 150 samples containing measurements of three distinct iris species. 
Ronald Fisher originally collected the majority of these samples in 1936 from the beautiful Gaspé Peninsula in Quebec!

This exercise was part of Udemy's Data Science and Machine Learning Bootcamp series provided by Pierian Data.

## Tuning of model hyperparameters
Here, the model was improved by performing a gridsearch to optimize hyperparameters. Note that training and cross-validation was performed on a distinct data set from testing, so the high score achieved by the optimized model should not indicate over-fitting on the training set. Indeed, from the initial data visualization it was clear to see that the three species were very distinct, which is why model prediction was very successful (and likely why this is a famous data set).
