clc
clear
close all
format long;
format compact;

fprintf('---------------------------------------------------------\n');
fprintf('      E-WOA: Enhanced Whale Optimization Algorithm\n')
fprintf('---------------------------------------------------------\n');


%% Problem Definition
% ----Modifiables----
bus_system = 33;            % [5, 33]
numDGs = 3;                 % [1, 2, 3] Recommended Maximum is 3             
numHours = 96;              % [1, 3, 24, 96]
maxDGsize_input = [];       % [DG1_Max, DG2_Max_, DG3_Max, DGN_Max]
Population = 50;
minDGSize = 0;  
MaxIt =  200;
run = 100;
%--------------------
   
minLocation = 2;            % Bus 1 is slack bus
VALID_HOURS = [1, 3, 24, 96];

if bus_system == 5
    maxLocation = 5;
    numDGs = 1;            
    numHours = 1;
end

if bus_system == 33
    maxLocation = 33;
end

if ~ismember(numHours, VALID_HOURS)
   error('Unsupported number of hours. Please choose 1, 3, 24, or 96.');
end

numVars = numVarsCalc(numHours, numDGs);
[BDactive, BDreactive, maxDGSize_per_hour] = loadBusData(numHours,bus_system);

if ~isempty(maxDGsize_input)
    maxDGSize = maxDGsize_input;
else
    maxDGSize = maxDGSize_per_hour;
end

% Variable Upper and Lower Bounds
upperbound = zeros(1, numDGs + numHours*numDGs);
upperbound(1:numDGs) = repmat(maxLocation, 1, numDGs);

for hour = 1:numHours
    for dg = 1:numDGs
        currentMaxSize = maxDGSize_per_hour(hour);
        if ~isempty(maxDGsize_input) && dg <= length(maxDGsize_input)
            if maxDGsize_input(dg) > currentMaxSize
                upperbound(numDGs + (hour-1)*numDGs + dg) = currentMaxSize;
            else
                upperbound(numDGs + (hour-1)*numDGs + dg) = maxDGsize_input(dg);
            end
        else
            upperbound(numDGs + (hour-1)*numDGs + dg) = currentMaxSize;
        end
    end
end

lowerbound = [repelem(minLocation, numDGs) repelem(minDGSize, numHours*numDGs)];


% Continue from previous optimization progress or start afresh

baseFilename = sprintf('EWOA_Results_%dBus_Hour%d_Test.mat', bus_system, numHours);
resultVarName = sprintf('EWOA_Result_%dBus_Hour%d', bus_system, numHours);

if exist(baseFilename, 'file')
    loadedData = load(baseFilename);
    completedRuns = size(loadedData.(resultVarName).Fitness, 2);
else
    completedRuns = 0;
  
    if bus_system == 33
        switch numHours
            case {1, 3, 24, 96}
                % For dynamic naming based on bus_system and numHours
                resultVarName = sprintf('EWOA_Result_%dBus_Hour%d', bus_system, numHours);
                optimalValueVarName = sprintf('EWOA_Optimal_value_%dBus_Hour%d', bus_system, numHours);
                optimalPositionVarName = sprintf('EWOA_Optimal_value_position_%dBus_Hour%d', bus_system, numHours);
                loadFlowResultsVarName = sprintf('EWOA_Load_Flow_Results_%dBus_Hour%d', bus_system, numHours);
                columnIndexVarName = sprintf('EWOA_column_index_%dBus_Hour%d', bus_system, numHours);

                loadedData.(resultVarName) = [];
                loadedData.(optimalValueVarName) = [];
                loadedData.(optimalPositionVarName) = [];
                loadedData.(loadFlowResultsVarName) = [];
                loadedData.(columnIndexVarName) = [];
        end
    end
end

fprintf('Progress: %3d%%\n',0)

% Main Loop
for j = (completedRuns + 1): run
    
    Runtime_Start = tic;
    Results = EWOA(Population,MaxIt,numVars,numDGs,numHours,maxDGSize,lowerbound,upperbound,BDactive,BDreactive);
    Runtime_Run = toc(Runtime_Start);
    
    loadedData.(resultVarName).Convergence(j,:) = Results.Convergence;
    loadedData.(resultVarName).Fitness(j)      = Results.Best_Fitness;
    loadedData.(resultVarName).Position(j,:)   = Results.Best_Position;
    loadedData.(resultVarName).Runtime(j,1)    = Runtime_Run;
    loadedData.(resultVarName).Runtime(j,2)    = Results.Runtime_Objective_Function;

    % Save the updated data after each run
    save(baseFilename, '-struct', 'loadedData');
    fprintf('\b\b\b\b\b%3.0f%%\n', j/run*100);
end


fprintf('\n');

optimalValueVarName = sprintf('EWOA_Optimal_value_%dBus_Hour%d', bus_system, numHours);
optimalPositionVarName = sprintf('EWOA_Optimal_value_position_%dBus_Hour%d', bus_system, numHours);
loadFlowResultsVarName = sprintf('EWOA_Load_Flow_Results_%dBus_Hour%d', bus_system, numHours);
columnIndexVarName = sprintf('EWOA_column_index_%dBus_Hour%d', bus_system, numHours);

%-----------Modify for Test System---------
% linedata33rds;
linedata5rds;
%------------------------------------------


if completedRuns == 100
    [Optimal_value_EWOA, column_index_EWOA] = min((loadedData.(resultVarName).Fitness));
    Optimal_value_position_EWOA = loadedData.(resultVarName).Position(column_index_EWOA,:);
    Load_Flow_Results_EWOA = Load_Flow_Result(Optimal_value_position_EWOA,numDGs,numHours,BDactive,BDreactive,linedata);

    loadedData.(optimalValueVarName) = Optimal_value_EWOA;
    loadedData.(optimalPositionVarName) = Optimal_value_position_EWOA;
    loadedData.(loadFlowResultsVarName) = Load_Flow_Results_EWOA;
    loadedData.(columnIndexVarName) = column_index_EWOA;
    save(baseFilename, '-struct', 'loadedData');
end


%% --------------Local Functions--------------

function formatted_num = three_decimal_format(n) 
    str = num2str(n,'%.3f');
    formatted_num = strrep(str, 'e+0', 'e');
end

function displayEWOAResults(EWOA_Result, run, Optimal_value_position_EWOA, numDGs, numHours)

    fprintf('---------------------------------------------------------\n');
    fprintf('The results are summarized as follows: \n');
    fprintf('    Number of runs = %d\n', run);
    fprintf('    Average = %s\n', three_decimal_format(mean(EWOA_Result.Fitness)));
    fprintf('    Min (Optimal Solution Value) = %s kW\n', three_decimal_format(min(EWOA_Result.Fitness)));
    fprintf('---------------------------------\n');

    locations = Optimal_value_position_EWOA(1:numDGs);
    dispatches = reshape(Optimal_value_position_EWOA(numDGs+1:end), numDGs, numHours);

    % Calculate the maximum dispatch for each DG
    maxDispatches = max(dispatches, [], 2);

    fprintf('DG\tLocation\tSize\n');
    fprintf('---------------------------------\n');
    for i = 1:numDGs
        fprintf('%d\t%d\t\t%.3f kW\n', i, locations(i), maxDispatches(i));
    end
    
    fprintf('---------------------------------\n');
end


function numVars = numVarsCalc(numHours, numDGs)
    numVars = (numHours * numDGs) + numDGs;
end

function [BDactive, BDreactive, maxDGSize] = loadBusData(numHours,bus_system)
    switch numHours
        case {1}
            if bus_system == 5
                busdata5radial;
            end
            if bus_system == 33
                load('peakactivedemand.mat', "BDactive");
                load('peakreactivedemand.mat', "BDreactive");
            end
        case {3}
            load('lowmidpeakactivedemand.mat', "BDactive");
            load('lowmidpeakreactivedemand.mat', "BDreactive");
        case {24}
            load('dailyactivedemand.mat', "BDactive");
            load('dailyreactivedemand.mat', "BDreactive");
        case {96}
            load('seasonalactivedemand.mat', "BDactive");
            load('seasonalreactivedemand.mat', "BDreactive");
    end
    maxDGSize = sum(BDactive);
end
