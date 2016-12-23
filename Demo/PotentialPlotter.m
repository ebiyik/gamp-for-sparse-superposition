clearvars; close all;

%% parameters
% general
B = 2; % section length
R = 0.3; % design rate
channel = 'bsc'; % awgnc | bec | bsc | zc
noiseParam = 0.1; % snr for awgnc, epsilon for bec, bsc or zc
firstE = 2^-5; % the first E value (MSE value) for which the potential will be calculated
lastE = 2^-0; % a little bit more than the last E value (MSE value) for which the potential will be calculated (not exact because potential function fails for exact 1)
exponentialStepSizeE = 0.01; % at each iteration, E value will be multiplied by 2^exponentialStepSizeE
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
