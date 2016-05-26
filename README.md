# Burghardt Models
Seed dormancy and seed dispersal models from Burghardt et al. 2015

To use the models create an instance, adjust any parameters that differ from the default and then use the run method.

For example:
```
% get the environmental data - at hourly intervals
envData   = csvread('environmental_data_Valencia_60yrs.csv',1,1);
moisture  = envData(:,5);
temp      = envData(:,4);

% create instance of the model
model = BurghardtGerminationModel();

% set some parameters to values different from the default
model.nSeedClasses = 21;

% run the model with the environment data and the starting index into that data
[germinationDay,fractionOfSeeds] = model.run(temp,moisture,1);
```
