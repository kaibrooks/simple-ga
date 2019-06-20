% Simple Genetic Algorithm
% Kai Brooks
% github.com/kaibrooks
% 2019
%
% A simple genetic algorithm example using MATLAB.
%
% Problem: find a solution to maximize x^2 for an n-bit number
%
% This is a trivial problem with a direct solution we can solve directly
% (not an NP problem). This program is meant to illustrate fundamentally
% how a GA finds solutions, using a simple problem for illustration.
%
% This code is not optimized, in the sense that many lines could be
% combined into a single function and the code minimized. However, since
% this is a basic algorithm, many lines are independent for clarity.


% Init
clc
close all
clear all
format
rng('default')


%% Options

popSize = 32;               % number of attempts per generation
maxGenerations = 2000;      % algorithm will run this many times
chromLength = 256;          % length of each chromosome, represents greatest binary value
f = @(x) x.^2;              % function to maximize


fmax = (bi2de(ones(1,chromLength),'left-msb'))^2; % find theoretical max (since we know it)


parentRate = 0.5;           % portion of members designated to breed each generation
eliteRate = 0.1;            % portion of members who are 'elites', and always breed (as a portion of total population, not parents)
mutationRate = 0.01;        % portion of members who mutate each generation

maxLoops = popSize*maxGenerations*chromLength;             % if any loop runs this many times with no changes, abort program. this prevents stalling

mutationNum = ceil(popSize*mutationRate); % number of members to mutate

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
    [~, eliteIndex] = maxk(fitness,eliteNum);  % find highest fitnesses, ~ removes saving the value since we don't need it
    
    isParent(eliteIndex) = 1;                  % mark it in isParent
    
    disp = [isParent fitness c];
    fitnessHistory;
    
    % -- Select the rest via roulette
    % Normalize fitnesses from 0 to 1
    normFitness = fitness - min(fitness(:));
    normFitness = normFitness ./ max(normFitness(:));
    
    normFitness(isnan(normFitness))=0.999;
    normFitness = normFitness + 0.001; % add a small number in case a bunch of pops have zero fitness
    
    % Subtract each members fitness from a random number, when < 0, select that
    % member (this ensures higher fitnesses have higher probilities).
    % This method starts at member 1 instead of picking a random position
    % because the order of members is already random
    
    miss = 0;  % check number of times loop doesn't do anything
    missMain = 0;
    
    while sum(isParent) <= parentNum    % pick more members until parentNum reached
        
        a = rand();             % random number from 0 to 1
        
        for j = 1:popSize
            
            if isParent(j), continue, end   % skip already selected
            
            a = a - normFitness(j);          % subtract members fitness
            
            if a < 0
                isParent(j) = 1;
                break;
            end
            
            if j == popSize     % back to start if no member found. this happens when a = rand() is larger than the sum of the remaining fitnesses
                j = 1;
                miss = miss+1;
            end
            
            if miss > maxLoops % kill program if loop stalls
                error('Hit maxLoops during roulette selection loop - for j=1:popSize')
            end
            
            
        end % roulette parent selection
        
        missMain = missMain+1;
        if missMain > maxLoops % kill program if loop stalls
            error('Hit maxLoops during roulette selection loop - while sum(isParent)')
        end
        
    end % roulette population selection
    
    [isParent fitness c];
    
    % Select two parents
    % Breed: Create new pop
    
    parentIndex = find(isParent);    % create vector of parents indexes
    newc = zeros(1,chromLength);
    
    % select two parents to breed
    for i = 1:2:parentNum
        parenta = c(parentIndex(i),:);
        parentb = c(parentIndex(i+1),:);
        
        for k = 1:4 %TODO fix number - each parent produces two offspring
            % select crossover point
            crossover = (randi(chromLength)-1); % pick a random point on the chromosome
            
            for j=1:chromLength
                if j < crossover
                    child(j) = parenta(j);  % first part
                else
                    child(j) = parentb(j);  % second part
                end
            end % end chromisome splicing
            newc = [newc; child];
            
            
        end % end two offspring each
        
    end
    
    % delete the temp entry
    newc(1, :) = [];    % cleanup the temp entry
    
    % mutate
    % bit flip mutation
    for i = 1:mutationNum
        randChromIndex = randi(chromLength);
        randPopIndex = randi(popSize);
        
        if newc(randPopIndex,randChromIndex) == 0
            newc(randPopIndex,randChromIndex) = 1;
            %fprintf('Mutated member %i at position %i - 0->1\n',randChromIndex,randPopIndex)
        else
           newc(randPopIndex,randChromIndex) = 0;
           %fprintf('Mutated member %i at position %i - 1->0\n',randChromIndex,randPopIndex)
        end
        
        
        
    end
    
        
    % write new chromosomes
    c = newc(randperm(size(newc,1)),:); % randomize the population (to ensure random roulette selection)
    
    
    fprintf('Generation: %4i \t Fitness: %e \t Optimized: %6.2f%%  \n',generation,fitnessHistory(generation), 100*(fitnessHistory(generation)/fmax))
end

plot(fitnessHistory)
xlabel('Generation')
ylabel('Fitness')
ylim([fmax*0.8 fmax])
grid on
% Shuffle

% Results: Plot
% Results: Display





% Functions
% Func: Binary to Decimal and back
% end 1:generation loop
