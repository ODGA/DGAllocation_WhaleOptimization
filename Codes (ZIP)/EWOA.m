% =========================================================================
%     Enhanced whale optimization algorithm (E-WOA) source codes 
%  Authors: Mohammad H.Nadimi-Shahraki, Hoda Zamani,Seyedali Mirjalili
%           --------------------------------------------- 
% Papers: Enhanced whale optimization algorithm for medical feature selection: A COVID-19 case study
%         Computers in Biology and Medicine,Volume 148, September 2022, 105858.
% https://www.sciencedirect.com/science/article/pii/S0010482522006126 
% https://doi.org/10.1016/j.compbiomed.2022.105858
%           ---------------------------------------------
% More materials are available on ResearchGate and for any questions please contact the authors.
% Emails:nadimi.mh@gmail.com, zamanie.hoda@gmail.com, ali.mirjalili@gmail.com
% =========================================================================
% This source code has been modified for use in solving the optimal
% allocation of distributed generation. 
% 
% Modifications:
%   1. Initialization: 
%           - Initialization of discrete variables
%           - Vectorized initialization of continuous variables
%   2. Added variable for measuring the runtime in processing the objective functions
%   3. Added local function to enforce boundary constraints
%   4. Checked boundaries and converted continuous variables to discrete 
%   variables after pool creation
%   5. Converted updated positions before fitness evaluation through
%   relaxation
%   6. Modified output variables for efficient storage
% =========================================================================

function Results = EWOA(Pop,MaxIter,numVars,numDGs,numHours,maxDGSize,lowerbound,upperbound,BDactive,BDreactive)

fobj = @Obj_Function; 
Convergence = zeros(1,MaxIter);

t = 0;

% --------Additional Data for Runtime-----------
Results.Runtime_Objective_Function = 0;

Positions = zeros(Pop, numVars);

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

Fitness = zeros(Pop, 1);
Objective_Runtime_Start = tic; 

for j = 1:Pop
    Fitness(j) = fobj(Positions(j,:), numDGs, numHours, BDactive, BDreactive);
end

Total_Runtime_Objective_Function = toc(Objective_Runtime_Start); 


% --------Initialization of Leader and Best Position-----------
[temp_fit, sorted_index] = sort(Fitness, 'ascend');
X_best = Positions(sorted_index(1),:);
Score_best = temp_fit(1);

% ---------------Initialization of Population Pool--------------- %%
P_rate =  20;                                   % The portion rate
Pool.Kappa = 1.5 * Pop;                         % The maximum pool size  
Pool.position = zeros(0, numVars);              % The solutions stored in the pool

worstPositionStartIndex = round(Pop - (Pool.Kappa * 0.3)+1);
X_worst = Positions(sorted_index(worstPositionStartIndex: end),:);
Pool = Pooling_Mechanism(Pool,numVars,X_worst,X_best);
% ------Boundary Check and Relaxation of Discrete Variables
for i = 1: size(Pool.position,1)
    Pool.position(i,:)   = simpleBoundConstraint(Pool.position(i,:), lowerbound,upperbound);
end
Pool.position(:,1:numDGs)   = round(Pool.position(:,1:numDGs),0);

t = 1;

% Main loop
while t <= MaxIter
    
    a2 = -1+t*((-1)/MaxIter);
   
    % Randomly select a P portion of the humpback whales
    P_portion = randperm(Pop,P_rate);
    
    Pop_Pool = Pool.position;
    [P_rnd1, P_rnd2] = RandIndex(size(Pop_Pool,1),  Pop);
    
    %% Probability rate
    p = rand(Pop, 1);
    
    %% Cauchy distribution
    A = 0.5 + 0.1 * tan(pi * (rand(Pop, 1) - 0.5));
    Idx = find(A <= 0);
    while ~ isempty(Idx)
        A(Idx) = 0.5 + 0.1 * tan(pi * (rand(length(Idx), 1) - 0.5));
        Idx = find(A <= 0);
    end
    A = min(A, 1);   
    
    for i = 1:size(Positions,1)
        if (i ~= P_portion)
            C = 2*rand;      % Eq. (6) in the paper
            b = 1;            
            for j=1:size(Positions,2)
                l=(a2-1)*rand+1;                
                if p(i) < 0.5
                    if A(i) < 0.5     % Enriched encircling prey search strategy
                        rand_leader_index = floor(size(Pop_Pool,1)*rand()+1); % Randomly selected from the matrix Pool
                        P_rnd3 = Pop_Pool(rand_leader_index, j);
                        D_prim = abs(C * X_best(j) - P_rnd3);  % Eq. (17)
                        X(i,j) =  X_best(j)- A(i)*D_prim;      % Eq. (16)
                        
                    elseif A(i)>= 0.5 % Preferential selecting search strategy
                        X(i,j)   = Positions(i,j) + A(i)*( C* Pop_Pool(P_rnd1(i),j)- Pop_Pool(P_rnd2(i),j)) ; % Eq. (15)
                    end
                elseif p(i)>=0.5      % Spiral bubble-net attacking
                    D_prim =  abs(X_best(j)- Positions(i,j));           % Eq. (10)
                    X(i,j)=  D_prim * exp(b.*l).*cos(2*pi*l)+ X_best(j);% Eq. (9)
                end
            end
        end
    end
    %% Migrating search strategy
    best_max = max(X_best);
    best_min = min(X_best);
    P_portion = P_portion';
    X_rnd  =  rand(size(P_portion,1),numVars).* (upperbound - lowerbound) + lowerbound;  % Eq.(13)
    X_brnd =  rand(size(P_portion,1),numVars).*(best_max - best_min)+ best_min  ; % Eq.(14)
    X(P_portion, :)  = X_rnd -  X_brnd; % Eq.(12)
    
    %% Computing the fitness of whales
    for i = 1: size(Positions,1)
        X(i,:)   = boundConstraint(X(i,:), Positions(i,:), lowerbound,upperbound);
        X(i,1:numDGs)   = round(X(i,1:numDGs),0);
        
        Objective_Runtime_Start = tic;
        Fit(i,1) = fobj(X(i,:),numDGs,numHours,BDactive,BDreactive);
        Runtime_Objective_Function = toc(Objective_Runtime_Start);
        Total_Runtime_Objective_Function = Total_Runtime_Objective_Function + Runtime_Objective_Function;
    

        if Fit(i,1) < Score_best
            Score_best =  Fit(i,1);
            X_best   =  X(i,:);
        end
    end    
    I = (Fitness > Fit);
    Ind = find(I == 1);
    X_worst = Positions(Ind, :);
    % Definition 1. (Pooling mechanism)
    Pool = Pooling_Mechanism(Pool, numVars, X_worst, X_best);
    for i = 1: size(Pool.position,1)
        Pool.position(i,:)   = simpleBoundConstraint(Pool.position(i,:), lowerbound,upperbound);
    end
    Pool.position(:,1:numDGs)   = round(Pool.position(:,1:numDGs),0);
    % Updating the position of whales
    Fitness(Ind) = Fit(Ind) ;
    Positions(Ind, :) = X(Ind, :);    
    
    Convergence(t) = Score_best;
    t = t + 1;
end

Results.Runtime_Objective_Function = Total_Runtime_Objective_Function;
Results.Convergence = Convergence;
Results.Best_Fitness = Score_best;        
Results.Best_Position = X_best;    

end

% -----------------Population Pool----------------------
function Pool = Pooling_Mechanism(Pool,problem_size,X_worst,X_best)
 
% Definition 1. (Pooling mechanism) 
best_max = max(X_best);
best_min = min(X_best);

B = randi([0 1], size(X_worst,1),size(X_worst,2));
% B = randi([0, 1], size(X_worst));
X_brnd = rand(size(X_worst,1),problem_size).*(best_max - best_min) + best_min ; % Eq.(14)
P = B.* X_brnd  + ~B .* X_worst;   % Eq.(11)


%% Store to Pool
PoolCrntSize = size(Pool.position,1);
FreeSpace = Pool.Kappa - PoolCrntSize ;
P = unique(P, 'rows');
if FreeSpace ~= 0
    if size (P,1) <= FreeSpace
        PopAll = [Pool.position ; P];
        PopAll = unique(PopAll, 'rows');
        Pool.position = PopAll;
    else
        Pool.position (PoolCrntSize+1:Pool.Kappa,:) = P(1:FreeSpace,:);
        Rm  = size(P,1) - FreeSpace;
        rnd = randi(PoolCrntSize,Rm,1);
        Pool.position(rnd,:) = P(FreeSpace+1:end,:);
    end
else
    rnd = randi(PoolCrntSize,size(P,1),1);
    Pool.position(rnd,:) = P;
end
Pool.position = unique(Pool.position, 'rows');
end

% -----------------Random Index----------------------
function [rnd1,rnd2] = RandIndex(NP1,SearchAgents_no)
warning off

rnd0 = 1 : SearchAgents_no;
NP0 = length(rnd0);

% == 1
rnd1  = floor(rand(1, NP0) * NP1) + 1;
Idx = (rnd1 == rnd0);
while sum(Idx) ~= 0
    rnd1(Idx) = floor(rand(1, sum(Idx)) * NP1) + 1;
    Idx = (rnd1 == rnd0);
end
% == 2
rnd2  = floor(rand(1, NP0) * NP1) + 1;
Idx = ((rnd2 == rnd1) | (rnd2 == rnd0));
while sum(Idx) ~= 0
    rnd2(Idx) = floor(rand(1, sum(Idx)) * NP1) + 1;
    Idx = ((rnd2 == rnd1) | (rnd2 == rnd0));
end
end

% ---------------Boundary Constraint-------------------
function X = boundConstraint (X, Positions, lb,ub)

% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound
%
% Version: 1.1   Date: 11/20/2007
% Written by Jingqiao Zhang, jingqiao@gmail.com
lu(1, :) = lb; lu(2, :) = ub;
[NP, ~] = size(Positions);   

%% check the lower bound
X_l = repmat(lu(1, :), NP, 1);
pos = X < X_l;
X(pos) = (Positions(pos) + X_l(pos)) / 2;

%% check the upper bound
xu = repmat(lu(2, :), NP, 1);
pos = X > xu;
X(pos) = (Positions(pos) + xu(pos)) / 2;
end

% -------Boundary Constraint for Pool---------- 
function X = simpleBoundConstraint(X, lb, ub)
X = max(X, lb);
X = min(X, ub);
end