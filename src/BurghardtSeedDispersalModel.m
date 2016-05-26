classdef BurghardtSeedDispersalModel < handle
    %BurghardtSeedDispersalModel Seed dispersal model from Burghardt 2015
    %   Simple seed dispersal model from Burghardt et al 2015, "Modeling
    %   the influence of genetic and environmental variation on the
    %   expression of plant life cycles across landscapes", The American
    %   Naturalist, 2015.
    %   
    %   This is a simple thremal unit summation model with a threshold
    %   for seed dispersal.
        
    properties
        threshold = 8448;  % parameter: thermal units threshold for disperal (default 8448)
        T_b       = 3;     % parameter: base temperature in degrees Celsuis (default 3)  
        
        storeProgress = false; % store hourly TU progress data (default false)
        progressData  = [];    % output hourly TU progess data (if stored, default [])
        
    end % properties
    
    methods
        function obj =  BurghardtSeedDispersalModel()
            % Constructs a model with default parameters. 
        end
        
        function dispersalTime = run(obj, temperature, start)
            % Runs the model and returns the dispersal time as the index
            % to the temperature data at which the threshold was first 
            % reached.
            % 
            % temperature: hourly temperature data in degrees Celsius
            % start: index into temperature array corresponding to flowering time
            %
            % returns the index into the temperature array at which
            % dispersal occurs. Or -1 if dispersal does not occur within
            % the given data range.
            
            obj.progressData = [];
            
            total = 0;
            
            for i = start:length(temperature)
                if (temperature(i) > obj.T_b)
                    total = total + (temperature(i)-obj.T_b);
                end
                
                if obj.storeProgress
                    obj.progressData = [obj.progressData;total];
                end
                    
                if (total >= obj.threshold) 
                    dispersalTime = i;
                    return;
                end
            end
            
            % did not reach the dispersal threshold
            dispersalTime = -1; 
            
        end % run function
        
    end % methods
    
end % classdef

