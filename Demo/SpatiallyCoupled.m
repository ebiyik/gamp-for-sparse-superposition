clearvars; close all;

%% parameters
% general
B = 32; % section length
R = 2.2; % design rate
L = 1024; % the number of sections
channel = 'awgnc'; % awgnc | bec | bsc | zc
noiseParam = 100; % snr for awgn, epsilon for bec, bsc or zc
n = 50; % the number of iterations 
% spatially coupling
betaSeed = 1.4; % the ratio between the height of seed and another block (only when the seed is larger first block)
Lr = 9; % the number of row blocks
Lc = 8; % the number of column blocks (must divide B*L)
wLeft = 7; % the number of window blocks at left
wRight = 1; % the number of window blocks at right
Jright = 0.2; % variance of the window blocks are right. (Jleft is 1 by default)

% for decoder
useHadamard = true; % true for using Hadamard-based operators, false otherwise
decoderInstances = 3; % the results will be averaged over this-many runs, 0 to disable decoder (it may be problematic to have large number of instances for BEC, BSC and ZC, as they occasionally fail due to NaNs, Infs, etc.)
useMse = false; % true to use average MSE for evaluation, false to use SER
saveMemory = false; % true for generating only the result, false for all iterations. It will be set to false if decoderShowPlot is true.

% result saving
saveResults = false; % true for saving, false for not saving as a mat-file
showPlot = true; % true for showing, false for no plot

%%%%% A PRIORI NO NEED TO CHANGE ANYTHING BELOW THIS LINE %%%%%

%% settings
rng('shuffle');
addpath(genpath('../Codes'));
showPlot = decoderInstances > 0 && showPlot;
saveMemory = ~showPlot && saveMemory;
if useMse; measure = 'Average MSE';
else measure = 'SER'; end;

%% Seed and bulk rates
N = B*L;
M = round(L*log2(B)/R);
subsectionRowSizes(1) = floor(L*log2(B)*betaSeed/(Lc*R));
subsectionRowSizes(1) = M-round((M-subsectionRowSizes(1))/(Lr-1))*(Lr-1);
subsectionRowSizes(2:Lr) = repmat((M-subsectionRowSizes(1))/(Lr-1),[Lr-1 1]);
Rsections(1) = L*log2(B)/(Lc*subsectionRowSizes(1));
Rsections(2:Lr) = L*log2(B)/(Lc*subsectionRowSizes(2));
disp(['Seed Rate:' num2str(Rsections(1))]);
disp(['Bulk Rate:' num2str(Rsections(2))]);

%% Decoder
fprintf('Decoder runs are starting...\n');
x = sparseSuperposition(B,L); % creates a random SS code
decoderResults = 0;
blockResults = 0;

for instanceNo=1:decoderInstances % Results will be averaged over instances
    fprintf('\tRun %d...', instanceNo);
    %% encoding
    if useHadamard
        J = generateSpatiallyCouplingMatrix(wLeft, wRight, Lr, Lc, Jright);
        rp = createRandomLinesAndSignsPermutationForOperators(Lc, Lr, J, subsectionRowSizes, N/Lc); 
        y = MultSeededHadamard(x', J, Lr, Lc, subsectionRowSizes, N/Lc, rp, []);

        %% set power P = 1 for each section
        rescale = zeros(Lr,1);
        for sectionNo=1:Lr
            rescale(sectionNo) = sqrt(mean(y(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo))).^2));
            y(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo))) = y(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo)))/rescale(sectionNo);
            J(sectionNo,:) = J(sectionNo,:)/(rescale(sectionNo)^2);
        end
    else
        J = generateSpatiallyCouplingMatrix(wLeft, wRight, Lr, Lc, Jright);
        % generate coding matrix A by using variance matrix J
        Acell = cell(Lr,Lc);
        for rowNo=1:Lr
            for colNo=1:Lc
                Acell{rowNo,colNo} = randn(subsectionRowSizes(rowNo),N/Lc)*sqrt(J(rowNo,colNo)/L);
            end
        end
        A = cell2mat(Acell);
        y = A * x;

        %% set power P = 1 for each section
        rescale = zeros(Lr,1);
        for sectionNo=1:Lr
            rescale(sectionNo) = sqrt(mean(y(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo))).^2));
            y(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo))) = y(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo)))/rescale(sectionNo);
            A(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo)),:) = A(sum(subsectionRowSizes(1:sectionNo-1))+1:sum(subsectionRowSizes(1:sectionNo)),:)/rescale(sectionNo);
            J(sectionNo,:) = J(sectionNo,:)/(rescale(sectionNo)^2);
        end
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
    blockTemp = zeros(size(xhat,2),Lc);
    if useMse
        for j=1:size(xhat,2)
            temp(j) = mse(x, xhat(:,j))*B;
        end
        for j=1:size(xhat,2)
            for k=1:Lc
                blockTemp(j,k) = mse(x((k-1)*N/Lc+1:k*N/Lc),xhat((k-1)*N/Lc+1:k*N/Lc,j))*B;
            end
        end
    else
        for j=1:size(xhat,2)
            temp(j) = ser(x, xhat_hardDecision(:,j), B);
        end
        for j=1:size(xhat,2)
            for k=1:Lc
                blockTemp(j,k) = ser(x((k-1)*N/Lc+1:k*N/Lc),xhat((k-1)*N/Lc+1:k*N/Lc,j), B);
            end
        end
    end
    blockResults = blockResults + blockTemp;
    decoderResults = decoderResults + temp; % will be averaged over instances

    fprintf(' Done!\n');
end
decoderResults = decoderResults / decoderInstances;
blockResults = blockResults / decoderInstances;

disp('Done!');

%% show/save results
if decoderInstances > 0
    disp(['Decoder Results - ' measure]);
    disp(decoderResults);
end

Legend = cell(0,0);
if showPlot
    semilogy(0:n, decoderResults(1:end), '-o');
    Legend{length(Legend)+1} = ['R = ' num2str(R), ' GAMP Decoder'];
    legend(Legend);
    title(['#Iterations vs ' measure ' for ' upper(channel) ' for B=' num2str(B) ', L=' num2str(L)]);
    xlabel('#Iterations');
    ylabel(measure);
    figure;
    for blockNo = 1:Lc
        semilogy(blockResults(:, blockNo), '-o');
        hold on;
        Legend{blockNo} = ['Block #' num2str(blockNo)];
    end
    legend(Legend);
    title('Block-by-Block Figure');
    xlabel('#Iterations');
    ylabel(measure);
end

if saveResults
    fileName = strcat('SC_', channel, num2str(noiseParam), '_B', num2str(B), '_R', num2str(R(1)), '_L', num2str(L), '_', num2str(decoderInstances), 'inst.mat');
    save(['results/' fileName]);
end
