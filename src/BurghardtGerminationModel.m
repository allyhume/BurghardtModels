classdef BurghardtGerminationModel < handle
    %BurghardtGerminationModel Seed germination model from Burghardt 2015
    %   Seed germination model from Burghardt et al 2015, "Modeling
    %   the influence of genetic and environmental variation on the
    %   expression of plant life cycles across landscapes", The American
    %   Naturalist, 2015.
    %   
        
    properties
        threshold = 1000;  % parameter: threshold for germination (default 1000)
        
        % Temperature
        
        T_bg      =    3;  % parameter: base temperature for germination in degrees Celsuis (default 3)
        T_o       =   22;  % parameter: optimal temperature for germination in degress Celsuis (default 22)
        k_T       =  0.12; % parameter: dormancy increase for each degree Celsuis above T_o (default 0.12)
        
        % Initial dormancy
        
        psi_mean       =  0; % parameter: mean dormancy at dispersal (default 0)
        psi_min        = -1; % parameter: mimumum dormancy possible (default -1)
        psi_breadth    =  1; % parameter: difference between lowest and highest dormancy classes (default 1)
        nSeedClasses   = 11; % parameter: number of seed dormancy classes (default 11)
        
        % Afterripening
        
        T_bar     = 3;      % parameter: base temperature for afterripening (default 3)
        psi_max   = -5;     % parameter: maximum moisture for afterripening (default -5)
        psi_l     = -350;   % parameter: lower moisture limit for afterripening (default -350)
        psi_u     = -50;    % parameter: upper moisture limit for afterripening (default -50)
        d_sat     = 40;     % parameter: days from 0 psi0_b to -1 psi_b (default 40)
        psi_scale = 1;
        
        % storing intermediate results
        
        storePsi = false; % store hourly psi data (default false)
        storeHTU = false; % store hourly HTU progress data (default false)
        
        psiData  = [];    % output hourly psi data if stored (default [])
        htuData  = [];    % output hourly HTU data if stored (default [])
        
    end % properties
    
    methods
        
        function obj =  BurghardtGerminationModel()
            % Constructs a model with default parameters. 
        end
        
        function [germinationDays,fractionOfSeeds] = run(obj, temperature, moisture, start)
            % Runs the model and returns arrays of germination days and
            % the number of fraction of seeds that germinate on that day.
            % Days are counted in 24 hour periods from the specified start
            % hour. The first day is day 1.
            % 
            % temperature: hourly temperature data in degrees Celsius
            % moisture: hourly temperature data in MPa
            % start: index into temperature and moisture arrays
            %        to start the model.
            %
            % returns arrays for the germination day and the fraction
            % of seeds that germinate on that day. If not all seeds 
            % germinate within the given environmental data time range
            % then the total fraction will not sum to 1.
            
            obj.psiData = [];
            obj.htuData = [];
            
            % set up initial psi values over the specified range
            psi = linspace(obj.psi_mean-(obj.psi_breadth/2), obj.psi_mean+(obj.psi_breadth/2), obj.nSeedClasses);
            
            % HTU and germination day for each seed class
            htu            = zeros(1,obj.nSeedClasses);
            germinationDay = zeros(1,obj.nSeedClasses);
    
            % Fraction of seeds for each seed class 
            nSeeds = normpdf(linspace(-3,3,obj.nSeedClasses),0,1);
            nSeeds = nSeeds/sum(nSeeds);  % total seeds is 1
            
            for t = start:length(temperature)
                [psi,htu] = obj.nextGerm(temperature(t),moisture(t),psi,htu);
        
                % set germination day for any that reach the threshold
                condition = (germinationDay == 0 & htu > obj.threshold);
                germinationDay = germinationDay + condition .* (fix(t/24) + 1);
                
                % store htu or psi values if required
                if obj.storeHTU 
                    obj.htuData = [obj.htuData;htu];
                end
                if obj.storePsi
                    obj.psiData = [obj.psiData;psi];
                end
                
                % break out if all have germinated
                if all(germinationDay > 0)
                    break;
                end
                
            end % t (time) for loop
            
            % Calc the fraction of seeds that germinate on each day
            germinationDays = unique(germinationDay);
            fractionOfSeeds = zeros(1,length(germinationDays));
            for i = 1:length(germinationDays)
                fractionOfSeeds(i) = sum((germinationDay==germinationDays(i)).*nSeeds);
            end
            
        end % run function
        
    end % public methods
    
    methods (Access = private)
        
        function [newPsi, newHTU] = nextGerm(obj,temperature,moisture,psi,htu)
            % Performs the next germination time step
            % 
            % temperature: temperature value at timestep
            % moisture: moisture value at timestep
            % psi: array of current psi values for the seed classes
            % htu: array of current HTU values for the seed classes
            %
            % returns the new psi and HTU values after the timestep
            
            newPsi = psi;
            newHTU = htu;
            
            isWet = (moisture > obj.psi_max);

            if (not(isWet)) 
                newPsi = obj.arLoss(temperature, moisture,psi);
            else
                newHTU = newHTU + obj.calcHTU(temperature, moisture, psi);
            end
            
        end % nextGerm function


        function newPsi = arLoss(obj,temperature, moisture, psi)
            % Calculates new psi values using the afterripening model
            %
            % temperature: temperature at current timestep
            % moisture: moisture at current timestep
            % psi: array of current psi values for the seed classes
            %
            % returns the new psi values after the timestep
            
            arSaturate = BurghardtGerminationModel.arHTU(20,-200,3,-350,-50)*24*obj.d_sat;
            newPsi = psi-(BurghardtGerminationModel.arHTU(temperature,moisture,obj.T_bar,obj.psi_l,obj.psi_u)*(obj.psi_scale/arSaturate));    
            newPsi = max(newPsi,obj.psi_min);
        end % arLoss function

    
        function addedHTU = calcHTU(obj, temperature, moisture, psi)
            % Calculates added HTU for this timestep
            %
            % temperature: temperature at current timestep
            % moisture: moisture at current timestep
            % psi: array of current psi values for the seed classes
            %
            % returns the HTU values for this timestep

            % Sub-optimal temperatures
            condition = psi<moisture & obj.T_bg < temperature & temperature <= obj.T_o;
            addedHTU = ((moisture-psi) * (temperature-obj.T_bg)) .* condition;

            %Supra-optimal temperatures
            mPsi = psi + obj.k_T*(temperature-obj.T_o);
            condition = temperature > obj.T_o & mPsi < moisture;
            addedHTU = addedHTU + (((moisture-mPsi) * (obj.T_o-obj.T_bg)) .* condition);

        end % calcHTU
        
    end % private methods
    
    methods (Static, Access = private)

        function x = arHTU(temperature, moisture, T_bar, psi_l, psi_u)
            % Calculates the afterripening HTU for one timestep
            %
            % temperature: temperature at this timestep
            % moisture: moisture at this timestep
            % T_bar: base temperature for afterripening
            % psi_l: lower moisture limit for afterripening 
            % psi_u: upper moisture limit for afterripening 
            %
            % returns the afterrippening HTU for the timestep
            x = 0;
            if moisture>=psi_l && moisture<=psi_u && temperature > T_bar
                x = ((psi_l-moisture)/(psi_l-psi_u))*(temperature-T_bar);
            end
            
            if moisture > psi_u && temperature > T_bar
                x = temperature-T_bar;
            end
        end % arHTU function
        
    end % static methods

end % classdef



