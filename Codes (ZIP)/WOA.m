%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes demo 1.0               %                                                             %
%  Developed in MATLAB R2011b(7.13)                                       %                                                       %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %                                                          %
%       Homepage: http://www.alimirjalili.com                             %                                                             %
%   Main paper: S. Mirjalili, A. Lewis                                    %
%               The Whale Optimization Algorithm,                         %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008   %                                                            %
%_________________________________________________________________________%
% This source code has been modified for use in solving the optimal       %
% allocation of distributed generation.                                   %
%                                                                         %
% Modifications:                                                          %
%   1. Initialization:                                                    %
%           - Initialization of discrete variables                        %
%           - Vectorized initialization of continuous variables           %
%   2. Added variable for measuring the runtime in processing the         %
% objective functions                                                     %
%   3. Converted updated positions before fitness evaluation through      %
%   relaxation                                                            %
%   4. Modified output variables for efficient storage                              %
% ________________________________________________________________________%

function Results = WOA(Pop,MaxIter,numVars,numDGs,numHours,maxDGSize,lowerbound,upperbound,BDactive,BDreactive)


fobj = @Obj_Function; 
Convergence = zeros(1,MaxIter);
Positions = zeros(Pop, numVars);
t = 0;

% --------Additional Data for Runtime-----------
Results.Runtime_Objective_Function = 0;
Total_Runtime_Objective_Function = 0;

% --------Leader Position and Score Initialization--------
Leader_pos = zeros(1,numVars);
Leader_score = inf;               

%---------Integer variables initialization-----------
for i = 1:numDGs
    Positions(:,i) = randi([lowerbound(i), upperbound(i)], Pop, 1);
end

%---------Continuous variables initialization------------
if numVars > numDGs
    Positions(:,numDGs+1:end) = bsxfun(@plus, lowerbound(numDGs+1:end), ...
                                       bsxfun(@times, rand(Pop, numVars-numDGs), ...
                                              (upperbound(numDGs+1:end) - lowerbound(numDGs+1:end))));
end


%% -----------WOA Main Loop-------------
while t < MaxIter
    
    % -----------Boundary Check-----------
    beyondUpperBound = Positions > upperbound;
    beyondLowerBound = Positions < lowerbound;
    Positions = (Positions .* ~(beyondUpperBound | beyondLowerBound)) + upperbound .* beyondUpperBound + lowerbound .* beyondLowerBound;
    Positions(:,1:numDGs) = round(Positions(:,1:numDGs),0);

    for i = 1:Pop
        Objective_Runtime_Start = tic;
        % --------Fitness Evaluation--------
        fitness = fobj(Positions(i,:),numDGs,numHours,BDactive,BDreactive);
        Runtime_Objective_Function = toc(Objective_Runtime_Start);
        Total_Runtime_Objective_Function = Total_Runtime_Objective_Function + Runtime_Objective_Function;
          
        % ----------Leader Update-----------
        if fitness < Leader_score             
            Leader_score = fitness;          
            Leader_pos = Positions(i,:);
        end
    end
    
    % -------Linearly Decreasing Parameters--------
    a = 2-t*((2)/MaxIter);              % Decreases from 2 to 0
    a2 = -1+t*((-1)/MaxIter);           % Decreases from -1 to -2
    
    % --------Position Update--------
    for i = 1:size(Positions,1)
        
        % --------Internal Parameters Update--------
        r1 = rand();          % Random number in [0,1]
        r2 = rand();          % Random number in [0,1]
        p = rand();           % Random number in [0,1]
        A = 2*a*r1-a;         % Encircling mechanism coefficent vector
        C = 2*r2;             % Encircling mechanism coefficent vector
        b = 1;                % Spiral model parameters
        l = (a2-1)*rand+1;    % Spiral model parameters
        
        for j = 1:size(Positions,2)
            if p < 0.5   
                if abs(A) >= 1
                    % ----------Search for Prey Strategy----------
                    rand_leader_index = floor(Pop*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand = abs(C*X_rand(j)-Positions(i,j));       
                    Positions(i,j) = X_rand(j) - A*D_X_rand;            
                    
                elseif abs(A) < 1
                    % ----------Encircling Prey Strategy----------
                    D_Leader = abs(C*Leader_pos(j) - Positions(i,j)); 
                    Positions(i,j) = Leader_pos(j) - A*D_Leader;        
                end
                
            elseif p >= 0.5
                % ----------Bubble-Net Hunting Strategy----------
                distance2Leader = abs(Leader_pos(j) - Positions(i,j));
                Positions(i,j) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + Leader_pos(j);
            end
        end
    end
    t = t + 1;
    Convergence(t) = Leader_score;
end

% ------------Output------------
Results.Runtime_Objective_Function = Total_Runtime_Objective_Function;
Results.Convergence = Convergence;
Results.Best_Fitness = Leader_score;        
Results.Best_Position = Leader_pos;    

end