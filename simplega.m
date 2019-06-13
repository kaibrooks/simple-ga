% Simple Genetic Algorithm
% Kai Brooks
% github.com/kaibrooks
% 2019
%
% A simple genetic algorithm example using MATLAB.
%
% Problem: find a solution to maximize x^2 for an 8-bit number
%
% This is a trivial problem with a direct solution we can solve directly
% (not an NP problem). This program is meant to illustrate fundamentally
% how a GA finds solutions, using a simple problem for illustration.
%
% This code is not optimized, in the sense that many lines could be
% combined into a single function and the code minimized. However, since
% this is a basic algorithm, each line is independent for clarity.


% Init
clc
close all
clear all
format
rng('shuffle')


%% Options

popSize = 8;                % number of attempts per generation
maxGenerations = 50;        % algorithm will run this many times
chromLength = 6;            % length of each chromosome, represents greatest binary value
f = @(x) x.^2;              % function to maximize

parentRate = 0.5;           % portion of members designated to breed each generation
eliteRate = 0.1;            % portion of members who are 'elites', and always breed (as a portion of total population, not parents)
mutationRate = 0.01;        % portion of members who mutate each generation

maxLoops = 100;             % if any loop runs this many times with no changes, abort program. this prevents stalling

% temp vars to wipe for release

%% Generate Initial Population
% ** Generate initial population
for i = 1:popSize
    c(:,i) = randi([0 1],chromLength,1);
end
c = transpose(c);   % flip it so each horizontal row represents a member

%% Solve and Compute Fitness (loop starts here)
for generation = 1:maxGenerations
    % Compute each member in decimal for input to the function
    decVals = bi2de(c,'left-msb'); % compute decimal values, leftmost MSB (...8 4 2 1)
    
    % Compute the function
    fitness = f(decVals);    % compute values using our defined function to get fitness scores
    fitnessHistory(generation) = mean(fitness);  % store each generations average fitness to see changes over time
    
    %% Breed New Generation
    % -- Breed: Select Parents
    % Get whole numbers for elites and parents
    parentNum = ceil(parentRate * popSize);     % total number of parents
    eliteNum = ceil(eliteRate * popSize);       % number of elites
    
    % If we have a non-even number of parents, add one to make it even
    if rem(parentNum, 2) ~= 0
        parentNum = parentNum + 1;
    end
    
    % Create a binary array to determine if a member will be selected to breed
    isParent = zeros(popSize,1);
    
    % Select elites first
    [~, selectedIndex] = maxk(fitness,eliteNum);  % find highest fitnesses, ~ removes saving the value since we don't need it
    isParent(selectedIndex) = 1;                  % mark it in isParent
    
    disp = [isParent fitness c];
    fitnessHistory
    
    % -- Select the rest via roulette
    % Normalize fitnesses from 0 to 1
    normFitness = fitness - min(fitness(:));
    normFitness = normFitness ./ max(normFitness(:));
    
    % Subtract each members fitness from a random number, when < 0, select that
    % member (this ensures higher fitnesses have higher probilities).
    % This method starts at member 1 instead of picking a random position
    % because the order of members is already random
    
    while sum(isParent) < (parentNum - eliteNum)+1    % pick more members until parentNum reached
        miss = 0;                       % check number of times loop doesn't do anything
        a = rand();                     % random number from 0 to 1
        for j = 1:popSize
            
            if isParent(j), continue, end   % skip already selected
            
            a = a - normFitness(j);          % subtract members fitness
            
            %fprintf('%i : %.3f - %.3f \n',j,a,normFitness(j))
            
            if a < 0
                isParent(j) = 1;
                break;
            end
            
            if j == popSize     % back to start if no member found. this happens when a = rand() is larger than the sum of the remaining fitnesses
                j = 1;
                miss = miss+1;
            end
            
            if miss > maxLoops % kill program if loop stalls
                error('Hit maxLoops during roulette selection')
            end
            
            
        end % roulette parent selection
        
    end % roulette population selection
    
    [isParent fitness c]
    
    
    % Breed: Select Parents
    % Breed: Selection process (tournament, roulette, elites, etc)
    
    % Breed: Create new pop
    
    % Attempt solution...
    
    
    % Results: Plot
    % Results: Display
    
    
    
    
    
    % Functions
    % Func: Binary to Decimal and back
    
    
end
