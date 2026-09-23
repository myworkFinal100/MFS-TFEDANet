% Just load the noisy EEG signal for two sessions
load('all_cont.mat');


% Butterworth Filtering (BPF of cut-off frequency 0.5 to 60 Hz)
Fs = 160;             
lowCut = 0.5;
highCut = 60;
order = 4;
[b,a] = butter(order, [lowCut highCut]/(Fs/2), 'bandpass');
all_cont_filt  = cell(size(all_cont));
for subj = 1:length(all_cont)
    eeg = all_cont{subj};      
    eeg_filt = zeros(size(eeg), 'single');
    for ch = 1:size(eeg,1)
        eeg_filt(ch,:) = single(filtfilt(b, a, double(eeg(ch,:))));
    end
    all_cont_filt{subj} = eeg_filt;
end


% Fourier mode decomposition
all_noisy= all_cont_filt;
num_cells = size(all_noisy, 1);
freq_bands = [0 4 8 13 30 60];
cont = applyFDM(all_noisy, num_cells, 64, Fs, freq_bands, 'dct');


% generate clean_Deta_Theta and cont_delta_theta 
numSubjects = numel(cont);
numChannels = 64;
cont_theta_delta = cell(numSubjects,1);
for subj = 1:numSubjects
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


% Extract Features no need for clean version of the signal
contFeatures = extractFeaturesER_TVS_P2P_HA_Skew(cont_theta_delta, cont);


% Detect contaminated CDTRs(using Otsu-based thresholding)
X = contFeatures;
mu = mean(X,1);
sigma = std(X,[],1);
sigma(sigma==0) = 1;
Xn = (X - mu) ./ sigma;
weights= [0.414341, 0.042751, 0.083180, 0.127385, 0.332343]; % these are the weights obtained using RF regression model
weights =weights(:);
MFS = Xn*weights;
threshold = zeros(numSubjects,1);
predictedLabel = false(size(MFS));
for subj = 1:numSubjects
    idx = (subj-1)*numChannels + (1:numChannels);
    MFSsub = MFS(idx);
    threshold(subj) = computeMFSThreshold1(MFSsub, "otsu");
    predictedLabel(idx) = detectOAChannels(MFSsub, threshold(subj));
end
fprintf('Detected %d / %d channel-instances\n', sum(predictedLabel), numSubjects*numChannels);


% Extract detected contaminated channels and corresponding clean channels
detectedContSignals  = {};
subjectIndex = [];
channelIndex = [];
signalLength = [];
k = 1;
for subj = 1:numSubjects
    for ch = 1:numChannels
        if predictedLabel(k)
            detectedContSignals{end+1,1}  = cont_theta_delta{subj}(ch,:);
            subjectIndex(end+1,1) = subj;
            channelIndex(end+1,1) = ch;
            signalLength(end+1,1) = size(cont_theta_delta{subj},2);
        end
        k = k + 1;
    end
end
fprintf('Detected contaminated channels = %d\n',numel(detectedContSignals));


% Segment detected noisy and clean channels into 4-second windows
segmentLength = 4*Fs;      
XT = {};
segmentInfo = [];
for i = 1:length(detectedContSignals)
    noisySig = detectedContSignals{i};
    L = floor(length(noisySig)/segmentLength)*segmentLength;
    noisySig = noisySig(1:L);
    numSegments = L/segmentLength;
    for s = 1:numSegments
        idx1 = (s-1)*segmentLength + 1;
        idx2 = s*segmentLength;
        XT{end+1,1} = noisySig(idx1:idx2);
        segmentInfo(end+1,:) = [...
            i,...                  
            subjectIndex(i),...   
            channelIndex(i),...    
            s];                    
    end
end
fprintf('Total training segments = %d\n',numel(XT));


% load trained model on real MMI dataset
load('Net_For_Real.mat');


% Make predictions
YPred = predict(net_RealMMI, XT, MiniBatchSize=16);
% Plotting
plot(YPred{400, 1})
hold on
plot(XT{400, 1})
legend('Predicted', 'Noisy')


% Reconstruction of the signal
numSignals = length(detectedContSignals);
reconstructedOutput = cell(numSignals,1);
for i = 1:numSignals
    idx = find(segmentInfo(:,1) == i);
    [~, order] = sort(segmentInfo(idx,4));
    idx = idx(order);
    recSig = [];
    for k = 1:length(idx)
        seg = YPred{idx(k)};
        seg = seg(:);
        recSig = [recSig; seg];
    end
    reconstructedOutput{i} = recSig';
end
fprintf('Reconstructed %d cleaned signals\n', numSignals);


% Place reconstructed signals back into the original subject-channel structure
reconstructedThetaDelta = cont_theta_delta;
for i = 1:numel(reconstructedOutput)
    subj = subjectIndex(i);
    ch   = channelIndex(i);
    N = size(reconstructedThetaDelta{subj}, 2);
    recSig = reconstructedOutput{i};
    recSig = recSig(:).';     
    L = min(length(recSig), N);
    reconstructedThetaDelta{subj}(ch,1:L) = recSig(1:L);
end
fprintf('Reconstructed signals placed back into original subject/channel positions.\n');
plot(cont_theta_delta{2, 1}(2,:))
hold on
plot(reconstructedThetaDelta{2, 1}(2,:))
legend('contaminated', 'predicted')


% Reconstruct cont structure using reconstructed theta+delta
cont_reconstructed = cont;   
for subj = 1:numSubjects
    for ch = 1:numChannels
        td = reconstructedThetaDelta{subj}(ch,:).';
        cont_reconstructed{subj}{ch}{4} = single(td);
        cont_reconstructed{subj}{ch}{5} = zeros(size(td), 'single');
    end
end
fprintf('cont_reconstructed created successfully.\n');


% Convert cont_reconstructed (cell structure) to reconstructed_ (64 x N matrix)
numSubjects = numel(cont_reconstructed);
numChannels = numel(cont_reconstructed{1});
reconstructed_ = cell(numSubjects,1);
for subj = 1:numSubjects
    N = length(cont_reconstructed{subj}{1}{1});
    temp = zeros(numChannels, N, 'single');
    for ch = 1:numChannels
        gamma = cont_reconstructed{subj}{ch}{1};
        beta  = cont_reconstructed{subj}{ch}{2};
        alpha = cont_reconstructed{subj}{ch}{3};
        theta = cont_reconstructed{subj}{ch}{4};
        delta = cont_reconstructed{subj}{ch}{5};
        temp(ch,:) = (gamma + beta + alpha + theta + delta).';
    end
    reconstructed_{subj} = temp;
end
fprintf('reconstructed_ created successfully.\\n');

numSubjects = numel(cont);
numChannels = numel(cont{1});
cont_ = cell(numSubjects,1);
for subj = 1:numSubjects
    N = length(cont{subj}{1}{1});
    temp = zeros(numChannels, N, 'single');
    for ch = 1:numChannels
        gamma = cont{subj}{ch}{1};
        beta  = cont{subj}{ch}{2};
        alpha = cont{subj}{ch}{3};
        theta = cont{subj}{ch}{4};
        delta = cont{subj}{ch}{5};
        temp(ch,:) = (gamma + beta + alpha + theta + delta).';
    end
    cont_{subj} = temp;
end
fprintf('reconstructed_ created successfully.\\n');
plot(reconstructed_{2, 1}(23,1:800))
hold on
plot(cont_{2, 1}(23,1:800))
legend('Predicted EEG', 'Contaminated EEG')
