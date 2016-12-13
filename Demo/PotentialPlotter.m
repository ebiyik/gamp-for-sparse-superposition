clearvars; close all;

%% parameters
% general
B = 4; % section length
R = 1.55; % design rate
channel = 'awgnc'; % awgnc | bec | bsc
noiseParam = 15; % snr for awgn, epsilon for bec or bsc
firstE = 2^-10; % the first E value (MSE value) for which the potential will be calculated
lastE = 1; % a little bit more than the last E value (MSE value) for which the potential will be calculated (not exact because potential function fails for exact 1)
exponentialStepSizeE = 0.05; % at each iteration, E value will be multiplied by 2^exponentialStepSizeE
monteCarloIterationCount = 20000; % The number of Monte Carlo iterations for the cases where B > 2

% result saving and showing
saveResults = true; % true for saving, false for not saving as a mat-file
showResults = true; % true for showing, false for no plot

%%%%% A PRIORI NO NEED TO CHANGE ANYTHING BELOW THIS LINE %%%%%

%% settings and initializations
rng('shuffle');
addpath(genpath('../Codes'));
eValues = 2.^(log2(firstE):exponentialStepSizeE:log2(lastE)-(1e-10));
FuValues = zeros(length(eValues),1);

%% potential calculations
for i=1:length(eValues)
    FuValues(i) = potential(B, R, channel, noiseParam, eValues(i), monteCarloIterationCount);
end

%% show/save results
if showResults
    semilogx(eValues, FuValues);
    Legend = ['R = ' num2str(R)];
    legend(Legend);
    title(['Potential for ' upper(channel) ' with noise parameter=' num2str(noiseParam) ' for B=' num2str(B)]);
    xlabel('E');
    ylabel('F_u(E)');
end
if saveResults
    filename = strcat('Potential_', channel, num2str(noiseParam), '_B', num2str(B), '_R', num2str(R), '.mat');
    save(['results/' filename], 'eValues', 'FuValues', 'R', 'B', 'channel', 'noiseParam');
end
