%% For SS (semisynthetic) Dataset
close all;
clear all;
clc;
cd('location where the files are stored') %location where the files (Pure_Data.mat, HEOG.mat, and VEOG.mat) are stored. These files can be
                                            % downloaded from the link provided in the manuscript for SS dataset.

% Load data
load('Pure_Data.mat')
numberofSessions = "put the number of sessions you have"; % Put the number of sessions you hvae in the dataset. As in the case of SS dataset used in the manuscript we have 54 sessions.
all_clean = cell(numberofSessions, 1);
for i = 1:numberofSessions
    resVarName = sprintf('sim%d_resampled', i);
    all_clean{i} = eval(resVarName);
end

load('HEOG.mat')
all_hoeg = cell(numberofSessions, 1);
for i = 1:numberofSessions
    resVarName = sprintf('heog_%d', i);
    all_hoeg{i} = eval(resVarName);
end

load('VEOG.mat')
all_voeg = cell(numberofSessions, 1);
for i = 1:numberofSessions
    resVarName = sprintf('veog_%d', i);
    all_voeg{i} = eval(resVarName);
end

clearvars -except all_voeg all_hoeg all_clean
%% Butterworth Filtering
% Parameters
Fs = "mention the sampling frequency";             % Change to your sampling frequency
lowCut = 0.5;                                       % Change to your low cut-off frequency
highCut = 60;                                       % Change to your high cut-off frequency
order = 4;
% Design Butterworth Bandpass Filter
[b,a] = butter(order, [lowCut highCut]/(Fs/2), 'bandpass');
% Initialize output cells
all_clean_filt = cell(size(all_clean));
all_hoeg_filt  = cell(size(all_hoeg));
all_voeg_filt  = cell(size(all_voeg));
% Filter EEG (all sessions)
for subj = 1:length(all_clean)
    eeg = all_clean{subj};     
    eeg_filt = zeros(size(eeg), 'single');
    for ch = 1:size(eeg,1)
        eeg_filt(ch,:) = single(filtfilt(b, a, double(eeg(ch,:))));
    end
    all_clean_filt{subj} = eeg_filt;
end
% Filter HEOG
for subj = 1:length(all_hoeg)
    hoeg = double(all_hoeg{subj}); 
    all_hoeg_filt{subj} = single(filtfilt(b, a, hoeg));
end
% Filter VEOG
for subj = 1:length(all_voeg)
    voeg = double(all_voeg{subj});     
    all_voeg_filt{subj} = single(filtfilt(b, a, voeg));
end
disp('Filtering completed.');
clearvars -except all_voeg_filt all_hoeg_filt all_clean_filt Fs

%% Parameters
numChannels = "Put the number of channels"; % here write the number of channels we have
contaminationPercent = 0.60; % 60% of the channels are to be artificially contaminated
numContaminated = round(contaminationPercent * numChannels);   
SNR_range = -12:-1; % Range for the SNR values
rng(1);
all_noisy = cell(size(all_clean_filt));
contaminated_channels = cell(size(all_clean_filt));
assigned_SNR = cell(size(all_clean_filt));
for subj = 1:length(all_clean_filt)
    eeg  = double(all_clean_filt{subj});      
    heog = double(all_hoeg_filt{subj});      
    veog = double(all_voeg_filt{subj});      
    ocular = heog + veog;
    noisyEEG = eeg;
    contam_idx = randperm(numChannels, numContaminated);
    snr_values = SNR_range(randi(numel(SNR_range), 1, numContaminated));
    for k = 1:numContaminated
        ch = contam_idx(k);
        targetSNR = snr_values(k);
        cleanSig = eeg(ch,:);
        artifact = ocular;
        Ps = mean(cleanSig.^2);
        Pa = mean(artifact.^2);
        scale = sqrt(Ps / (Pa * 10^(targetSNR/10)));
        artifact_scaled = scale * artifact;
        noisyEEG(ch,:) = eeg(ch,:) + artifact_scaled;
    end
    all_noisy{subj} = single(noisyEEG);
    contaminated_channels{subj} = contam_idx;
    assigned_SNR{subj} = snr_values;
end
disp('Random 60% channel contamination completed successfully.');


segmentLength = 4 * Fs;    
for subj = 1:length(all_clean_filt)
    N = size(all_clean_filt{subj},2);
    Nnew = floor(N/segmentLength) * segmentLength;
    all_clean_filt{subj} = all_clean_filt{subj}(:,1:Nnew);
    all_noisy{subj}      = all_noisy{subj}(:,1:Nnew);
end
%% Fourier Mode Decomposition
all_clean = all_clean_filt;
all_noisy = all_noisy;
num_cells = size(all_noisy, 1);
freq_bands = [0 4 8 13 30 60];
cont = applyFDM(all_noisy, num_cells, numChannels, Fs, freq_bands, 'dct');
clean = applyFDM(all_clean, num_cells, numChannels, Fs, freq_bands, 'dct');
% generate clean_Deta_Thata and cont_delta_theta (combined delta theta rhythms CDTRs)
numberofSessions = numel(clean);
clean_theta_delta = cell(numberofSessions,1);
for subj = 1:numberofSessions
    channelData = clean{subj};   
    N = length(channelData{1}{4});  
    temp = zeros(numChannels, N, 'single');
    for ch = 1:numChannels
        theta = channelData{ch}{4};   
        delta = channelData{ch}{5};   
        temp(ch,:) = (theta + delta).';
    end
    clean_theta_delta{subj} = temp;
end

numberofSessions = numel(cont);
cont_theta_delta = cell(numberofSessions,1);
for subj = 1:numberofSessions
    channelData = cont{subj};
    N = length(channelData{1}{4});
    temp = zeros(numChannels, N, 'single');
    for ch = 1:numChannels
        theta = channelData{ch}{4};
        delta = channelData{ch}{5};
        temp(ch,:) = (theta + delta).';
    end
    cont_theta_delta{subj} = temp;
end
%%  Extract features from CDTRs
contFeatures = extractFeaturesER_TVS_P2P_HA_Skew(cont_theta_delta, cont);

X = contFeatures;
mu = mean(X,1);
sigma = std(X,[],1);
sigma(sigma==0) = 1;
Xn = (X - mu) ./ sigma;
weights= [0.414341, 0.042751, 0.083180, 0.127385, 0.332343]; % These are the weights computed through the RF-based regression model kept common for all datasets.
weights =weights(:);
MFS = Xn*weights;
threshold = zeros(numberofSessions,1);
predictedLabel = false(size(MFS));
for subj = 1:numberofSessions
    idx = (subj-1)*numChannels + (1:numChannels);
    MFSsub = MFS(idx);
    threshold(subj) = computeMFSThreshold1(MFSsub, "otsu");
    predictedLabel(idx) = detectOAChannels(MFSsub, threshold(subj));
end
%% Extract detected OA-contaminated channels and corresponding clean channels
detectedContSignals  = {};
detectedCleanSignals = {};
subjectIndex = [];
channelIndex = [];
signalLength = [];
k = 1;
for subj = 1:numberofSessions
    for ch = 1:numChannels
        if predictedLabel(k)
            detectedContSignals{end+1,1}  = cont_theta_delta{subj}(ch,:);
            detectedCleanSignals{end+1,1} = clean_theta_delta{subj}(ch,:);
            subjectIndex(end+1,1) = subj;
            channelIndex(end+1,1) = ch;
            signalLength(end+1,1) = size(cont_theta_delta{subj},2);
        end
        k = k + 1;
    end
end
fprintf('Detected contaminated channels = %d\n',numel(detectedContSignals));
%% Segment detected noisy and clean channels into 4-second windows
segmentLength = 4*Fs; 
XT = {};
YT = {};
% [SignalID SubjectID ChannelID SegmentNo]
segmentInfo = [];
for i = 1:length(detectedContSignals)
    noisySig = detectedContSignals{i};
    cleanSig = detectedCleanSignals{i};
    L = floor(length(noisySig)/segmentLength)*segmentLength;
    noisySig = noisySig(1:L);
    cleanSig = cleanSig(1:L);
    numSegments = L/segmentLength;
    for s = 1:numSegments
        idx1 = (s-1)*segmentLength + 1;
        idx2 = s*segmentLength;
        XT{end+1,1} = noisySig(idx1:idx2);
        YT{end+1,1} = cleanSig(idx1:idx2);
        segmentInfo(end+1,:) = [...
            i,...                  % Detected signal ID
            subjectIndex(i),...    % Subject
            channelIndex(i),...    % Channel
            s];                    % Segment number
    end
end
fprintf('Total training segments = %d\n',numel(XT));
%%  TFEDANet Model Training
XTn = cell(size(XT));
YTn = cell(size(YT));
muXT = zeros(numel(XT),1,'single');
sigmaXT = zeros(numel(XT),1,'single');
for k = 1:numel(XT)
    x = single(XT{k});
    y = single(YT{k});
    mu = mean(x);
    sigma = std(x);
    if sigma < 1e-8
        sigma = 1;
    end
    XTn{k} = (x - mu) / sigma;
    YTn{k} = (y - mu) / sigma;

    muXT(k) = mu;
    sigmaXT(k) = sigma;
end

%10-Fold Cross Validation
cv = cvpartition(numel(XTn),'KFold',10);
denoisedSegments = cell(size(XTn));
for fold = 1:10
    fprintf('\n=============================\n');
    fprintf('Fold %d / 10\n',fold);
    fprintf('=============================\n');
    % Split
    trainIdx = training(cv,fold);
    testIdx  = test(cv,fold);
    XTrain = XTn(trainIdx);
    YTrain = YTn(trainIdx);
    XTest = XTn(testIdx);
    XTrain = cellfun(@single,XTrain,'UniformOutput',false);
    YTrain = cellfun(@single,YTrain,'UniformOutput',false);
    XTest = cellfun(@single,XTest,'UniformOutput',false);
    % Validation split (10%)
    cvVal = cvpartition(numel(XTrain),'HoldOut',0.10);
    idxTrain = training(cvVal);
    idxVal   = test(cvVal);
    XVal = XTrain(idxVal);
    YVal = YTrain(idxVal);
    XTrain = XTrain(idxTrain);
    YTrain = YTrain(idxTrain);

    
% Model Training
winLength = 44;
overlapLength = 16;
numFeatures=66;
win=rectwin(winLength);
input = sequenceInputLayer(1, MinLength=800, Name="input");
BN = batchNormalizationLayer(Name="BN");
stft = stftLayer(Window=win, OverlapLength=overlapLength, transform="realimag", Name="stft");
convA1 = convolution1dLayer(32, 16, Padding="same", Stride=4, Name="convA1");
reluA1 = geluLayer(Name="reluA1");
convA2 = convolution1dLayer(32, 16, Padding="same", Stride=4, Name="convA2");
reluA2 = geluLayer(Name="reluA2");
convB1 = convolution1dLayer(32, 16, Padding="same", Stride=4, Name="convB1");
reluB1 = geluLayer(Name="reluB1");
transConvB2 = transposedConv1dLayer(32, 16, Cropping="same", Stride=4, Name="transConvB2");
reluTransB2 = geluLayer(Name="reluTransB2");
concatSkip = concatenationLayer(1, 2, Name="concatSkip");
transConvC1 = transposedConv1dLayer(32, 16, Cropping="same", Stride=4, Name="transConvC1");
reluTransC1 = geluLayer(Name="reluTransC1");
transConvC2 = transposedConv1dLayer(32, 16, Cropping="same", Stride=4, Name="transConvC2");
reluTransC2 = geluLayer(Name="reluTransC2");
flattenC2 = flattenLayer(Name= 'ft')
concat = concatenationLayer(1, 2, Name="concat");
bilstm = bilstmLayer(200, Name="bilstmA");
drop= dropoutLayer(0.4);
FC= fullyConnectedLayer(46, Name="fcMerge");
ISTFT= istftLayer(Window=win, OverlapLength=overlapLength, Name="istft");
FC2= fullyConnectedLayer(1, Name="fcFinal");
Regression = regressionLayer(Name="regressionOutput");
lgraph = layerGraph();
lgraph = addLayers(lgraph, input);
lgraph = addLayers(lgraph, BN);
lgraph = addLayers(lgraph, stft);
lgraph = addLayers(lgraph, [convA1; reluA1; convA2; reluA2]);         
lgraph = addLayers(lgraph, [convB1; reluB1; transConvB2; reluTransB2]);
lgraph = addLayers(lgraph, concatSkip);                                
lgraph = addLayers(lgraph, [transConvC1; reluTransC1; transConvC2; reluTransC2; flattenC2; concat; bilstm; drop; FC;ISTFT;FC2;Regression]); % Branch C
lgraph = connectLayers(lgraph, "input", "BN");
lgraph = connectLayers(lgraph, "BN", "stft");
lgraph = connectLayers(lgraph, "stft", "convA1");
lgraph = connectLayers(lgraph, "reluA2", "convB1");
lgraph = connectLayers(lgraph, "reluA2", "concatSkip/in1");
lgraph = connectLayers(lgraph, "reluTransB2", "concatSkip/in2");
lgraph = connectLayers(lgraph, "concatSkip", "transConvC1");
flatten = flattenLayer(Name='FT');
norm0 = layerNormalizationLayer(Name="norm0");
attention = selfAttentionLayer(8,128,'AttentionMask','causal', Name = "ATT");
fc1_1 = fullyConnectedLayer(256, Name="fc1_1");
relu1 = geluLayer(Name="relu1");
fc1_2 = fullyConnectedLayer(46, Name="fc1_2");
add1 = additionLayer(2, Name="add1");
norm1 = layerNormalizationLayer(Name="norm1");
attention2 = selfAttentionLayer(16,64, Name='attention2');
fc2_1 = fullyConnectedLayer(256, Name="fc2_1");
relu2 = geluLayer(Name="relu2");
fc2_2 = fullyConnectedLayer(46, Name="fc2_2");
add2 = additionLayer(2, Name="add2");
norm2 = layerNormalizationLayer(Name="norm2");
attention3 = selfAttentionLayer(8,128, Name='attention3');
fc3_1 = fullyConnectedLayer(256, Name="fc3_1");
relu3 = geluLayer(Name="relu3");
fc3_2 = fullyConnectedLayer(46, Name="fc3_2");
add3 = additionLayer(2, Name="add3");
norm3 = layerNormalizationLayer(Name="norm3");
FC = fullyConnectedLayer(16, Name="FC3");
lgraph = addLayers(lgraph, flatten);
lgraph = addLayers(lgraph, norm0);
lgraph = addLayers(lgraph, attention);
lgraph = addLayers(lgraph, fc1_1);
lgraph = addLayers(lgraph, relu1);
lgraph = addLayers(lgraph, fc1_2);
lgraph = addLayers(lgraph, add1);
lgraph = addLayers(lgraph, norm1);

lgraph = addLayers(lgraph, attention2);
lgraph = addLayers(lgraph, fc2_1);
lgraph = addLayers(lgraph, relu2);
lgraph = addLayers(lgraph, fc2_2);
lgraph = addLayers(lgraph, add2);
lgraph = addLayers(lgraph, norm2);
lgraph = addLayers(lgraph, attention3);
lgraph = addLayers(lgraph, fc3_1);
lgraph = addLayers(lgraph, relu3);
lgraph = addLayers(lgraph, fc3_2);
lgraph = addLayers(lgraph, add3);
lgraph = addLayers(lgraph, norm3);
lgraph = addLayers(lgraph, FC);
lgraph = connectLayers(lgraph, "stft", "FT");
lgraph = connectLayers(lgraph, "FT", "norm0");
lgraph = connectLayers(lgraph, "norm0", "ATT");
lgraph = connectLayers(lgraph, "ATT", "fc1_1");
lgraph = connectLayers(lgraph, "fc1_1", "relu1");
lgraph = connectLayers(lgraph, "relu1", "fc1_2");
lgraph = connectLayers(lgraph, "fc1_2", "add1/in1");
lgraph = connectLayers(lgraph, "norm0", "add1/in2");
lgraph = connectLayers(lgraph, "add1", "norm1");
lgraph = connectLayers(lgraph, "norm1", "attention2");
lgraph = connectLayers(lgraph, "attention2", "fc2_1");
lgraph = connectLayers(lgraph, "fc2_1", "relu2");
lgraph = connectLayers(lgraph, "relu2", "fc2_2");
lgraph = connectLayers(lgraph, "fc2_2", "add2/in1");
lgraph = connectLayers(lgraph, "norm1", "add2/in2");
lgraph = connectLayers(lgraph, "add2", "norm2");
lgraph = connectLayers(lgraph, "norm2", "attention3");
lgraph = connectLayers(lgraph, "attention3", "fc3_1");
lgraph = connectLayers(lgraph, "fc3_1", "relu3");
lgraph = connectLayers(lgraph, "relu3", "fc3_2");
lgraph = connectLayers(lgraph, "fc3_2", "add3/in1");
lgraph = connectLayers(lgraph, "norm2", "add3/in2");
lgraph = connectLayers(lgraph, "add3", "norm3");
lgraph = connectLayers(lgraph, "norm3", "FC3");
lgraph = connectLayers(lgraph, "FC3", "concat/in2");
analyzeNetwork(lgraph);
layers = lgraph;

options = trainingOptions("adam", ...
    MaxEpochs=50, ...
    MiniBatchSize=16, ...
    InitialLearnRate=0.001, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.5, ...
    LearnRateDropPeriod=20, ...
    ValidationData={XVal,YVal},...
    ValidationFrequency=20, ...
    ValidationPatience=40, ...
    L2Regularization=0.001, ...
    GradientThreshold=1, ...
    Plots="training-progress", ...
    Shuffle="every-epoch", ...
    Verbose=false);

    % Train the model
    net = trainNetwork(XTrain,YTrain,layers,options);
    % Predict the output
    YPred = predict(net, XTest, MiniBatchSize=16);
    % Denormalize and store prediction
    testIndices = find(testIdx);
    for n = 1:numel(YPred)
        idx = testIndices(n);
        denoisedSegments{idx} = YPred{n} * sigmaXT(idx) + muXT(idx);
    end
end
disp('10-Fold completed.')
%% Signal Reconstruction;
numDetected = length(subjectIndex);
reconstructedOutput = cell(numDetected,1);
for i = 1:numDetected
    idx = find(segmentInfo(:,1)==i);    
    [~,ord] = sort(segmentInfo(idx,4));
    idx = idx(ord);
    reconstructedOutput{i} = [];
    for k = 1:length(idx)
        reconstructedOutput{i} = ...
            [reconstructedOutput{i}, denoisedSegments{idx(k)}];
    end
end

reconstructed_theta_delta = cont_theta_delta;
for i = 1:numDetected
    subj = subjectIndex(i);
    ch   = channelIndex(i);
    reconstructed_theta_delta{subj}(ch,:) = reconstructedOutput{i};
end

% Replace only reconstructed channels
cont_reconstructed =cont;
for i = 1:numDetected
    subj = subjectIndex(i);
    ch   = channelIndex(i);
    cont_reconstructed{subj}{ch}{4} = ...
        reconstructed_theta_delta{subj}(ch,:).';
    cont_reconstructed{subj}{ch}{5} = ...
        zeros(size(cont_reconstructed{subj}{ch}{4}),'single');
end

% Reconstruct full EEG by summing all 5 FDM rhythms
numberofSessions = numel(cont_reconstructed);
reconstructed_ = cell(numberofSessions,1);
for subj = 1:numberofSessions
    N = length(cont_reconstructed{subj}{1}{1});
    eeg = zeros(numChannels,N,'single');
    for ch = 1:numChannels
        temp = zeros(N,1,'single');
        for r = 1:5
            temp = temp + cont_reconstructed{subj}{ch}{r};
        end
        eeg(ch,:) = temp.';
    end
    reconstructed_{subj} = eeg;
end

numberofSessions = numel(clean);
clean_ = cell(numberofSessions,1);
for subj = 1:numberofSessions
    N = length(clean{subj}{1}{1});
    eeg = zeros(numChannels,N,'single');
    for ch = 1:numChannels
        temp = zeros(N,1,'single');
        for r = 1:5
            temp = temp + clean{subj}{ch}{r};
        end
        eeg(ch,:) = temp.';
    end
    clean_{subj} = eeg;
end

numberofSessions = numel(cont);
cont_ = cell(numberofSessions,1);
for subj = 1:numberofSessions
    N = length(cont{subj}{1}{1});
    eeg = zeros(numChannels,N,'single');
    for ch = 1:numChannels
        temp = zeros(N,1,'single');
        for r = 1:5
            temp = temp + cont{subj}{ch}{r};
        end
        eeg(ch,:) = temp.';
    end
    cont_{subj} = eeg;
end
plot(cont_{10, 1}(1,1:800))
hold on
plot(reconstructed_{10, 1}(1,1:800))
hold on
plot(clean_{10, 1}(1,1:800))
legend('contaminated', 'predicted', 'groundTruth')