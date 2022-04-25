%
%
% Project Title: Advacned Time-Invariant Multi-Objective Particle Swarm Optimization (AT-MOPSO)
%
% Main Implementation of AT-MOPSO
% 
% 

clc;
clear;
close all;

%% Problem Definition

% Number of Objective Functions
nObj = 5;

% v = norm_query_probs
% w = pop_block_peers
% x = numPeers
% y = peer_weights
% z = numBlocks

CostFunction = @(v,w,x,y,z) Cost_Prob_Storage_Occupancy(v,w,x,y,z);      % Cost Function

nVar = 3;             % Number of Decision Variables == this will be Y in our case === the number of peers

VarSize = [1 nVar];   % Size of Decision Variables Matrix

numBlocks = 200;     % the number of blocks of each peer

VarMin = 0;          % Lower Bound of Variables
VarMax = 1;          % Upper Bound of Variables

F0 = 0.95;            % Initial query frequency of all  blocks
a1 = 0.01;           % Exp Decay for Exp Freq. case
a2 = 0.004;          % Linear Decay for Linear Freq. case

[sum_query_probs_of_blocks, query_probs_of_blocks] = CalculateQueryProbability_FFixed(F0, numBlocks);    % Get the query frequencies of all blocks -- Fixed freq. case 
%[sum_query_probs_of_blocks, query_probs_of_blocks] = CalculateQueryProbability_FLinear(F0, numBlocks, a2);    % Get the query frequencies of all blocks -- Linear freq. case 
%[sum_query_probs_of_blocks, query_probs_of_blocks] = CalculateQueryProbability_FExp(F0, numBlocks, a1);    % Get the query frequencies of all blocks -- Exp. freq. case 

peer_weights = [0.3 0.2 0.5];       % weights of all peers for their

obj_funcs_constraints = [0.45; 0.45; 0.45; 3.0; 0.6];    % the constraints for the value of the objective functions == always numBlocks + 2

obj_funcs_weights = [0.16; 0.16; 0.16; 0.12; 0.4];


%% MOPSO Parameters

MaxIt = 200;           % Maximum Number of Iterations

nPop = 100;            % Population Size

nRep = 100;            % Repository Size

w = 0.5;              % Inertia Weight
w1 = 0.7;
w2 = 0.3;
wdamp = 0.99;         % Intertia Weight Damping Rate
c1 = 1;               % Personal Learning Coefficient
c2 = 2;               % Global Learning Coefficient

nGrid = 7;            % Number of Grids per Dimension
alpha = 0.1;          % Inflation Rate

beta = 2;             % Leader Selection Pressure
gamma = 2;            % Deletion Selection Pressure

mu = 0.1;             % Mutation Rate


%% NSGA3 Params

% Colect Parameters
nDivision = 10;
Zr = GenerateReferencePoints(nObj, nDivision);

params.nPop = nPop;
params.Zr = Zr;
params.nZr = size(Zr, 2);
params.zmin = [];
params.zmax = [];
params.smin = [];

pMutation = 0.5;       % Mutation Percentage
nMutation = round(pMutation*nPop);  % Number of Mutants

mu = 0.02;     % Mutation Rate

sigma = 0.1*(VarMax-VarMin); % Mutation Step Size


%% Initialization

empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
empty_particle.Blocks = [];
empty_particle.IsDominated = [];
empty_particle.GridIndex = [];
empty_particle.GridSubIndex = [];
empty_particle.Rank = [];
empty_particle.DominationSet = [];
empty_particle.DominatedCount = [];
empty_particle.NormalizedCost = [];
empty_particle.AssociatedRef = [];
empty_particle.DistanceToAssociatedRef = [];

%% 

% averages for cost functions
averages.f1 = [];
averages.f2 = [];
averages.f3 = [];
averages.f4 = [];
averages.f5 = [];


pop = repmat(empty_particle, nPop, 1);

for i = 1:nPop
    
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    
    pop(i).Velocity = zeros(VarSize);
    
    %pop(i).Cost = CostFunction(pop(i).Position);
    [pop(i).Cost, pop(i).Blocks] = CostFunction(sum_query_probs_of_blocks, pop(i).Position, nVar, peer_weights, numBlocks);
    
    
    % Update Personal Best
    pop(i).Best.Position = pop(i).Position;
    pop(i).Best.Cost = pop(i).Cost;
    pop(i).Best.Blocks = pop(i).Blocks;
    
end

% Determine Domination
pop = DetermineDomination(pop);

rep = pop(~[pop.IsDominated]);

Grid = CreateGrid(rep, nGrid, alpha);

for i = 1:numel(rep)
    rep(i) = FindGridIndex(rep(i), Grid);
end


%% MOPSO Main Loop

for it = 1:MaxIt
    
    %{
    use_rand = randi(100);
    use_mod = mod(use_rand,2);
    
    if use_mod == 0
        c = 2.0 + (4.0-2.0)*rand();
        phi = c + c;
        CF = (2/abs(2 - phi - sqrt(phi^2 - (4*phi))));
    else
        w_t = ((w1-w2)*((MaxIt-it)/MaxIt)) + w2;
        c1 = 1;               
        c2 = 2;
    end
    %}
    
    for i = 1:nPop
        
        leader = SelectLeader(rep, beta);
        
        %{
        if use_mod == 0
            pop(i).Velocity = CF.*(pop(i).Velocity ...
                +c*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
                +c*rand(VarSize).*(leader.Position-pop(i).Position));
        else 
            pop(i).Velocity = w_t*pop(i).Velocity ...
                +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
                +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        end
        %}
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
        
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        pop(i).Position = max(pop(i).Position, VarMin);
        pop(i).Position = min(pop(i).Position, VarMax);
        
        %pop(i).Cost = CostFunction(pop(i).Position);
        [pop(i).Cost, pop(i).Blocks] = CostFunction(sum_query_probs_of_blocks, pop(i).Position, nVar, peer_weights, numBlocks);
    end
        
        %{
        % Apply Mutation
        pm = (1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
            NewSol.Position = Mutate(pop(i).Position, pm, VarMin, VarMax);
            NewSol.Cost = CostFunction(NewSol.Position);
            if Dominates(NewSol, pop(i))
                pop(i).Position = NewSol.Position;
                pop(i).Cost = NewSol.Cost;

            elseif Dominates(pop(i), NewSol)
                % Do Nothing

            else
                if rand<0.5
                    pop(i).Position = NewSol.Position;
                    pop(i).Cost = NewSol.Cost;
                end
            end
        end
        
        if Dominates(pop(i), pop(i).Best)
            pop(i).Best.Position = pop(i).Position;
            pop(i).Best.Cost = pop(i).Cost;
            
        elseif Dominates(pop(i).Best, pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position = pop(i).Position;
                pop(i).Best.Cost = pop(i).Cost;
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep = [rep
         pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep = DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep = rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid = CreateGrid(rep, nGrid, alpha);

    % Update Grid Indices
    for i = 1:numel(rep)
        rep(i) = FindGridIndex(rep(i), Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra = numel(rep)-nRep;
        for e = 1:Extra
            rep = DeleteOneRepMemebr(rep, gamma);
        end
        
    end
    
    % Plot Costs
    figure(1);
    PlotCosts(pop, rep);
    pause(0.01);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    %}    
    
    % Mutation
    popm = repmat(empty_particle, nMutation, 1);
    for k = 1:nMutation

        i = randi([1 nPop]);
        p = pop(i);

        popm(k).Position = Mutate(p.Position, mu, sigma);
        popm(k).Velocity = p.Velocity;
        popm(k).Best.Position = p.Best.Position;
        popm(k).Best.Cost = p.Best.Cost;

        %popm(k).Cost = CostFunction(popm(k).Position);
        [popm(k).Cost, popm(k).Blocks] = CostFunction(sum_query_probs_of_blocks, popm(k).Position, nVar, peer_weights, numBlocks);
        
        if Dominates(popm(k), popm(k).Best)
            pop(k).Best.Position = pop(k).Position;
            pop(k).Best.Cost = pop(k).Cost;
        end

    end

    % Merge
    pop = [pop
           popm]; %#ok
    
    % Sort Population and Perform Selection
    [pop, F, params] = SortAndSelectPopulation(pop, params);
    
    % Store F1
    F1 = pop(F{1});

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);

    %{
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1);
    pause(0.01);
    %}
    
    % get the averages of cost functions for plotting
    best_iter_costs = [pop.Cost];
    
    averages.f1(it) = mean(best_iter_costs(1,:));
    averages.f2(it) = mean(best_iter_costs(2,:));
    averages.f3(it) = mean(best_iter_costs(3,:));
    averages.f4(it) = mean(best_iter_costs(4,:));
    averages.f5(it) = mean(best_iter_costs(5,:));
    
    % Damping Inertia Weight
    w = w*wdamp;

    
end

%% Resluts

%% Objective fucntions to plot

plotAverages.f1 = [];
plotAverages.f2 = [];
plotAverages.f3 = [];
plotAverages.f4 = [];
plotAverages.f5 = [];

for plot_f = 1:numel(averages.f1)
    if plot_f == 1
        plotAverages.f1(1) = averages.f1(1);
        plotAverages.f2(1) = averages.f2(1);
        plotAverages.f3(1) = averages.f3(1);
        plotAverages.f4(1) = averages.f4(1);
        plotAverages.f5(1) = averages.f5(1);
    end
    
    if mod(plot_f,10) == 0
        arr_index = (plot_f/10) + 1;
        plotAverages.f1(arr_index) = averages.f1(plot_f);
        plotAverages.f2(arr_index) = averages.f2(plot_f);
        plotAverages.f3(arr_index) = averages.f3(plot_f);
        plotAverages.f4(arr_index) = averages.f4(plot_f);
        plotAverages.f5(arr_index) = averages.f5(plot_f);
    end
end
%% Results plotted for the objective functions
x = 0:10:200;

figure(1)
plot(x,plotAverages.f1, '-+');
figure(2)
plot(x,plotAverages.f2, '-o');
figure(3)
plot(x,plotAverages.f3, '-x');
figure(4)
plot(x,plotAverages.f4, '-*');
figure(5)
plot(x,plotAverages.f5, '-s');


%% Filter the solution set to get only acceptable results based on constraints

% a matrix to hold the filtered values
filtered_rep = repmat(empty_particle, 1, 1);

for filter = 1:numel(pop)
    if all(pop(filter).Cost <= obj_funcs_constraints) && any(pop(filter).Cost < obj_funcs_constraints)
        if numel(filtered_rep(1).Cost) == 0
            filtered_rep(1) = pop(filter);
            
        else
            filtered_rep = [filtered_rep
                            pop(filter)];
        end
        
    end
end


%% Finding the most suitable solution from filtered set using weighted sum

best_solution = empty_particle;

%set weighted sum of the costs of the best solution to inf
best_solution.w_sum = inf;

for best = 1:numel(filtered_rep)
    product_peer = obj_funcs_weights.*filtered_rep(best).Cost;
    weighted_sum_peer = sum(product_peer);
    
    if weighted_sum_peer < best_solution.w_sum
        best_solution.Position = filtered_rep(best).Position;
        best_solution.Blocks = filtered_rep(best).Blocks;
        best_solution.Cost = filtered_rep(best).Cost;
        best_solution.w_sum = weighted_sum_peer; 
    end
end

best_solution

