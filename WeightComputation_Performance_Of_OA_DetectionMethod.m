close all;
clear all;
clc;

% First download the dataset and then proceed
cd('location where SS dataset is stored')

% Load raw EEG data for SS dataset
load('Pure_Data.mat')
all_clean = cell(54, 1);
for i = 1:54
    resVarName = sprintf('sim%d_resampled', i);
    all_clean{i} = eval(resVarName);
end
load('HEOG.mat')
all_hoeg = cell(54, 1);
for i = 1:54
    resVarName = sprintf('heog_%d', i);
    all_hoeg{i} = eval(resVarName);
end
load('VEOG.mat')
all_voeg = cell(54, 1);
for i = 1:54
    resVarName = sprintf('veog_%d', i);
    all_voeg{i} = eval(resVarName);
end
clearvars -except all_voeg all_hoeg all_clean


% Butterworth Filtering
Fs = 'sampling rate of the EEG signal';             
lowCut = 0.5;
highCut = 60;
order = 4;
[b,a] = butter(order, [lowCut highCut]/(Fs/2), 'bandpass');
all_clean_filt = cell(size(all_clean));
all_hoeg_filt  = cell(size(all_hoeg));
all_voeg_filt  = cell(size(all_voeg));

for subj = 1:length(all_clean)
    eeg = all_clean{subj};      
    eeg_filt = zeros(size(eeg), 'single');
    for ch = 1:size(eeg,1)
        eeg_filt(ch,:) = single(filtfilt(b, a, double(eeg(ch,:))));
    end
    all_clean_filt{subj} = eeg_filt;
end

for subj = 1:length(all_hoeg)
    hoeg = double(all_hoeg{subj});      
    all_hoeg_filt{subj} = single(filtfilt(b, a, hoeg));
end

for subj = 1:length(all_voeg)
    voeg = double(all_voeg{subj});      % 1 x N
    all_voeg_filt{subj} = single(filtfilt(b, a, voeg));
end

clearvars -except all_voeg_filt all_hoeg_filt all_clean_filt Fs

% Mixing EEG VEOG and HEOG data to generate OA-contaminated EEG
numChannels = 'Enter the number of channels';
contaminationPercent = 0.60;
numContaminated = round(contaminationPercent * numChannels);   % 11 channels
SNR_range = -12:-1;
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
clearvars -except  all_noisy all_clean_filt Fs contaminated_channels assigned_SNR


segmentLength = 4 * Fs;    
for subj = 1:length(all_clean_filt)
    N = size(all_clean_filt{subj},2);
    Nnew = floor(N/segmentLength) * segmentLength;
    all_clean_filt{subj} = all_clean_filt{subj}(:,1:Nnew);
    all_noisy{subj}      = all_noisy{subj}(:,1:Nnew);
end

% Fourier Mode Decomposition for Rhythem Separation
all_clean = all_clean_filt;
all_noisy = all_noisy;
num_cells = size(all_noisy, 1);
freq_bands = [0 4 8 13 30 60];
cont = applyFDM(all_noisy, num_cells, 19, Fs, freq_bands, 'dct');
clean = applyFDM(all_clean, num_cells, 19, Fs, freq_bands, 'dct');

% CDTRs
numSession = numel(clean);
clean_theta_delta = cell(numSession,1);
for subj = 1:numSession
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

numSession = numel(cont);
cont_theta_delta = cell(numSession,1);
for subj = 1:numSession
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


% Extract features from CDTRs
contFeatures = extractFeaturesER_TVS_P2P_HA_Skew(cont_theta_delta, cont);

% Labeling
y = zeros(numSession*numChannels,1);
idx = 1;
for subj = 1:numSession
    contIdx = contaminated_channels{subj};
    for ch = 1:numChannels
        if ismember(ch,contIdx)
            y(idx) = 1;
        end
        idx = idx + 1;
    end
end

X = contFeatures;
mu = mean(X,1);
sigma = std(X,[],1);
sigma(sigma==0) = 1;
Xn = (X - mu) ./ sigma;
[model, weights] = trainMFSRegression(Xn, y, 'RF');
weights = weights(:)
MFS = Xn*weights;

% Otsu-based Thresholding and Prediction
threshold = zeros(numSession,1);
predictedLabel = false(size(MFS));
for subj = 1:numSession
    idx = (subj-1)*numChannels + (1:numChannels);
    MFSsub = MFS(idx);
    ysub   = y(idx);
    threshold(subj) = computeMFSThreshold1(MFSsub, "otsu");
    predictedLabel(idx) = detectOAChannels(MFSsub, threshold(subj));
end
predicted_channels = cell(numSession,1);
idx = 1;
for subj = 1:numSession
    predicted_channels{subj} = find(predictedLabel(idx:idx+18));
    idx = idx + 19;
end
TP = sum((predictedLabel == 1) & (y == 1));
TN = sum((predictedLabel == 0) & (y == 0));
FP = sum((predictedLabel == 1) & (y == 0));
FN = sum((predictedLabel == 0) & (y == 1));

accuracy    = (TP + TN) / (TP + TN + FP + FN);
precision   = TP / (TP + FP + eps);
recall      = TP / (TP + FN + eps);     
specificity = TN / (TN + FP + eps);
F1score     = 2 * precision * recall / (precision + recall + eps);

fprintf('Accuracy    = %.2f %%\n', accuracy*100);
fprintf('Precision   = %.2f %%\n', precision*100);
fprintf('Recall      = %.2f %%\n', recall*100);
fprintf('Specificity = %.2f %%\n', specificity*100);
fprintf('F1-score    = %.2f %%\n', F1score*100);
