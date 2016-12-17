clearvars; close all;

%% parameters
% general
B = 4; % section length
R = 1.3; % design rate
channel = 'awgnc'; % awgnc | bec | bsc
noiseParam = 15; % snr for awgn, epsilon for bec or bsc
n = 20; % the number of iterations 

% for decoder
L = 1024; % the number of sections
useHadamard = false; % true for using Hadamard-based operators, false otherwise
decoderInstances = 3; % the results will be averaged over this-many runs, 0 to disable decoder (it may be problematic to have large number of instances for BEC and BSC, as they occasionally fail due to NaNs, Infs, etc.)
decoderShowPlot = true; % true for showing, false for no plot
useMse = true; % true to use average MSE for evaluation, false to use SER
saveMemory = false; % true for generating only the result, false for all iterations. It will be set to false if decoderShowPlot is true.

% for state evolution (B, R, channel, noiseParam, n, useMse are from the decoder parameters above.)
stateEvolutionInstances = 3; % the results will be averaged over this many runs, 0 to disable state evolution (it may be problematic to have large number of instances for BEC and BSC, as they occasionally fail due to NaNs, Infs, etc.)
stateEvolutionShowPlot = true; % true for showing, false for no plot
monteCarloIterationCount = 20000;  % The number of Monte Carlo iterations

% result saving
saveResults = true; % true for saving, false for not saving as a mat-file

%%%%% A PRIORI NO NEED TO CHANGE ANYTHING BELOW THIS LINE %%%%%

%% settings
rng('shuffle');
addpath(genpath('../Codes'));
decoderShowPlot = decoderInstances > 0 && decoderShowPlot;
saveMemory = ~decoderShowPlot && saveMemory;
stateEvolutionShowPlot = stateEvolutionInstances > 0 && stateEvolutionShowPlot;
if useMse; measure = 'Average MSE';
else measure = 'SER'; end;

%% Decoder
fprintf('Decoder runs are starting...\n');
x = sparseSuperposition(B,L); % creates a random SS code
M = round(L*log2(B)/R);
N = B*L;
decoderResults = 0;
for instanceNo=1:decoderInstances % Results will be averaged over instances
    fprintf('\tRun %d...', instanceNo);
    %% encoding
    if useHadamard
        rp = createRandomLinesAndSignsPermutationForOperators(1, 1, 1/L, M, B*L); 
        y = MultSeededHadamard(x', 1/L, 1, 1, M, B*L, rp, []);

        %% set power P = 1 for each section
        rescale = mean(y.^2);
        y = y ./ sqrt(rescale);
    else
        A = randn(M,B*L);
        y = A * x;

        %% set power P = 1 for each section
        rescale = mean(y.^2);
        y = y ./ sqrt(rescale);
        A = A ./ sqrt(rescale);
    end

    %% communication over the channel
    switch channel
        case 'awgnc'
            y = y + randn(size(y)) * sqrt(1/noiseParam);
        case 'bsc'
            multipliers = ones(length(y(:)),1);
            multipliers(1:round(noiseParam*length(y(:)))) = -1;
            multipliers = multipliers(randperm(length(y)));
            y = (2*double(y>0)-1).*multipliers;
        case 'bec'
            multipliers = ones(length(y(:)),1);
            multipliers(1:round(noiseParam*length(y(:)))) = 0;
            multipliers = multipliers(randperm(length(y)));
            y = (2*double(y>0)-1).*multipliers;
    end

    %% GAMP decoder step
    % run the appropriate algorithm
    if useHadamard
        J = 1/L;
        Lr = 1;
        Lc = 1;
        subsectionRowSizes = M;
        if saveMemory
            gamp_memorySave_hadamard;
        else
            gamp_hadamard;
        end
    else
        if saveMemory
            gamp_memorySave;
        else
            gamp;
        end
    end

    %% hard decision
    xhat_hardDecision = hardDecision(xhat, B);

    %% result evaluations
    temp = zeros(size(xhat,2),1);
    if useMse
        for j=1:size(xhat,2)
            temp(j) = mse(x, xhat(:,j))*B;
        end
    else
        for j=1:size(xhat,2)
            temp(j) = ser(x, xhat_hardDecision(:,j), B);
        end
    end
    decoderResults = decoderResults + temp; % will be averaged over instances

    fprintf(' Done!\n');
end
decoderResults = decoderResults / decoderInstances;

%% State evolution
fprintf('State Evolution runs are starting...\n');
stateEvolutionResults = 0;
for instanceNo=1:stateEvolutionInstances % Results will be averaged over instances
    fprintf('\tRun %d...', instanceNo);
    [seTemp_mse, seTemp_ser] = stateEvolution(B, R, noiseParam, n, channel, monteCarloIterationCount);
    if useMse; stateEvolutionResults = stateEvolutionResults + seTemp_mse;
    else stateEvolutionResults = stateEvolutionResults + seTemp_ser; end
    fprintf(' Done!\n');
end
stateEvolutionResults = stateEvolutionResults / stateEvolutionInstances;

disp('Done!');

%% show/save results
if decoderInstances > 0
    disp(['Decoder Results - ' measure]);
    disp(decoderResults);
end
if stateEvolutionInstances > 0
    disp(['State Evolution Results - ' measure]);
    disp(stateEvolutionResults);
end

Legend = cell(0,0);
if decoderShowPlot
    semilogy(0:n, decoderResults, '-o');
    hold on;
    Legend{length(Legend)+1} = ['R = ' num2str(R), ' GAMP Decoder'];
    legend(Legend);
    title(['#Iterations vs ' measure ' for ' upper(channel) ' for B=' num2str(B) ', L=' num2str(L)]);
    xlabel('#Iterations');
    ylabel(measure);
end
if stateEvolutionShowPlot
    semilogy(0:n, stateEvolutionResults, '-*');
    hold off;
    Legend{length(Legend)+1} = ['R = ' num2str(R), ' State Evolution'];
    legend(Legend);
    title(['#Iterations vs ' measure ' for ' upper(channel) ' for B=' num2str(B) ', L=' num2str(L)]);
    xlabel('#Iterations');
    ylabel(measure);
end

if saveResults
    fileName = strcat('NonSC_', channel, num2str(noiseParam), '_B', num2str(B), '_R', num2str(R(1)), '_L', num2str(L), '_', num2str(decoderInstances), 'inst.mat');
    save(['results/' fileName]);
end
