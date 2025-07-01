These data and codes implement part of the analysis conducted in the paper

Samad, I., Sutaria, D. & Shanker, K. (in review). Spatio-temporal drivers of elasmobranch catch are site and fishery specific: Implications for fisheries management and species conservation. Biodiversity and Conservation.

The code is structured in three parts which include
a) Tuning each model to obtain optimal hyperparameters for the data
b) Run each model and generate predictions (including ensembles) while obtaining model evaluation metrics and thresholds
c) Generate evaluation metrics for ensemble models and CBI metrics for all models.

Note that different parameters of the setup e.g., buffer distance, thinning radius can be changed for estimating changes in model performance.