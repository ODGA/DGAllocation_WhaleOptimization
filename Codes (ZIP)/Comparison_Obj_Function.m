function fitness = Comparison_Obj_Function(x,numDGs,numHours,BDactive,BDreactive)   

%% Particle Position Extraction
DGsize = x(numDGs+1:end);
DGloc = x(1:numDGs);

% -----DG Location Constraint-----
sameLocation = numel(DGloc) ~= numel(unique(DGloc));
Penalty = 500000;
if sameLocation
    fitness = Penalty;
    return;
end

%% Line Data Input
linedata33rds_comparison; 
LD = linedata;

%% Bus Data Modification

for hour = 1:numHours
    for i1 = 1:numDGs
        index = (hour - 1)*numDGs+i1;
        BDactive(DGloc(i1),hour) = BDactive(DGloc(i1),hour) - DGsize(index);
    end
end

%% Base Values (Apparent Power, Voltage, Impedance)
Sbase = 100;                    % MVA
Vbase = 12.66;                  % kV
Zbase = (Vbase^2)/Sbase;      
%% Per Unit Conversion of Line Impedance and Load
LD(:,4:5) = LD(:,4:5)/Zbase;
BDactive = BDactive/(1000*Sbase);
BDreactive = BDreactive/(1000*Sbase);      
%% Complex Load and Impedance in PU
Sload = complex(BDactive,BDreactive);        
Z = complex(LD(:,4),LD(:,5));             
%% Initialization
V = ones(size(BDactive,1),1);                   
Iline = zeros(size(LD,1),1);           
MaxError = 1;
Tolerance = 0.0000001;
Iter = 0;
Voltage = [];
Ploss = [];

for hours = 1:numHours                             
    while (MaxError>Tolerance)
        Iter = Iter + 1;                                                            % Counter
        V0 = V;                                                                     % V_(iter-1) for Solving Maximum Error
        % -----Backward Sweep (Solving for Currents)-----
        Iload = conj(Sload(:,hours)./V);                                                     % Load Currents
        for j = size(LD,1):-1:1                                                     % From End to Starting Bus
            row_temp = [];
            % col_temp = [];
            [row_temp,~] = find(LD(:,2:3)==LD(j,3));                         % Position of Buses Connected to Current Bus
            if size(row_temp,1)==1                                                  % If Only One Bus is Connected then Current Bus is End or Starting Bus 
                Iline(LD(j,1)) = Iload(LD(j,3));                                    % Load Current of Bus is Equal to Line Current
            else
                Iline(LD(j,1)) = Iload(LD(j,3))+sum(Iline(LD(row_temp,1)))-Iline(LD(j,1));    % The Additional Third Term Subtracts Line Current from Previous Iteration
            end
        end
        % -----Forward Sweep (Solving for Voltages)-----
        for j = 1:size(LD,1)
             V(LD(j,3)) = V(LD(j,2))-Iline(LD(j,1))*Z(j);
        end
        % -----Maximum Error----- 
        MaxError = max(abs(V - V0));                                                % V_current_iter - V_previous_iter
    end    

    Voltage = [Voltage abs(V)];
    Ploss = [Ploss real(Z.*(abs(Iline.^2)))*100000];

    V = ones(size(BDactive,1),1);
    MaxError = 1;
    Iter = 0;
end

TPloss = sum(Ploss,"all");                           % Total P Loss kW

%% Constraints
% -----Voltage Constraints-----
Vmax = 1.0;
Vmin = 0.9;
VisBelowVmin = any(Voltage < Vmin);
VisAboveVmax = any(Voltage > Vmax);

%% Fitness Value
if any([VisBelowVmin, VisAboveVmax])
    fitness = Penalty;
    return;
end

fitness = TPloss;

end     