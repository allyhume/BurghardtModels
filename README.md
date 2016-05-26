# Burghardt Models
Seed dormancy and seed dispersal models from Burghardt et al. 2015

Matlab implementations based on the original R code from: https://github.com/lianaburghardt/Arabimodel_BurghardtAmNat2015

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

% visualise the results
figure;
bar(germinationDay, fractionOfSeeds);
```

The models have flag properties that allows you to specify that you wish to store some of the intermediate results.  

For the germination model you can store the hourly psi and HTU data.  For example:
```
model = BurghardtGerminationModel();
model.storePsi = true;
model.storeHTU = true;

[germinationDay,fractionOfSeeds] = model.run(temp,moisture,1);

% see the intermediate data
figure;
hold on
for ii = 1:nSeedClasses
  plot(model.htuData(:,ii))
end
hold off

figure;
hold on
for ii = 1:nSeedClasses
  plot(model.psiData(:,ii))
end
hold off
```

For the seed disperal model you can store the thermal units progress data.  For example:

```
model = BurghardtSeedDispersalModel();
model.storeProgress = true;
start = 3000;
dispersalTime = model.run(temp,start);

fprintf('Dispersal time: %d\n', dispersalTime);

figure;
plot(start:dispersalTime,model.progressData);
```

