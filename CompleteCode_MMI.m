
all_clean = "clean the MMI sessions using EWT-SCA method";
all_cont = "Noisy MMI sessions to be cleaned";
%% Butterworth Filtering
Fs = 160;           
lowCut = 0.5;
highCut = 60;
order = 4;
[b,a] = butter(order, [lowCut highCut]/(Fs/2), 'bandpass');
all_clean_filt = cell(size(all_clean));
all_cont_filt  = cell(size(all_cont));
for subj = 1:length(all_clean)
    eeg = all_clean{subj};
    eeg_filt = zeros(size(eeg), 'single');
    for ch = 1:size(eeg,1)
        eeg_filt(ch,:) = single(filtfilt(b, a, double(eeg(ch,:))));
    end
    all_clean_filt{subj} = eeg_filt;
end
for subj = 1:length(all_cont)
    cont = double(all_cont{subj});     
    all_cont_filt{subj} = single(filtfilt(b, a, cont));
end

segmentLength = 4 * Fs;    
for subj = 1:length(all_clean_filt)
    N = size(all_clean_filt{subj},2);
    Nnew = floor(N/segmentLength) * segmentLength;
    all_clean{subj} = all_clean_filt{subj}(:,1:Nnew);
    all_noisy{subj}      = all_cont_filt{subj}(:,1:Nnew);
end
%% Fourier decomposition method (FDM)
all_clean = all_clean;
all_noisy= all_noisy';
num_cells = size(all_noisy, 1);
freq_bands = [0 4 8 13 30 60];
cont = applyFDM(all_noisy, num_cells, 64, Fs, freq_bands, 'dct');
clean = applyFDM(all_clean, num_cells, 64, Fs, freq_bands, 'dct');
% generate CDTRs 
numSubjects = numel(clean);
numChannels = 64;
clean_theta_delta = cell(numSubjects,1);
for subj = 1:numSubjects
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
%% MFS
contFeatures = extractFeaturesER_TVS_P2P_HA_Skew(cont_theta_delta, cont);
% Detect contaminated CDTR channels (using Otsu-based thresholding)
X = contFeatures;
mu = mean(X,1);
sigma = std(X,[],1);
sigma(sigma==0) = 1;
Xn = (X - mu) ./ sigma;
weights= [0.414341, 0.042751, 0.083180, 0.127385, 0.332343]; % weights computed using RF-based regression model using SS dataset
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
% Extract detected OA-contaminated channels and corresponding clean channels
detectedContSignals  = {};
detectedCleanSignals = {};
subjectIndex = [];
channelIndex = [];
signalLength = [];
k = 1;
for subj = 1:numSubjects
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

segmentLength = 4*Fs;    

XT = {};
YT = {};
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
            i,...                  
            subjectIndex(i),...   
            channelIndex(i),...   
            s];                   
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  TFEDANet (10-Fold Cross Validation) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cv = cvpartition(numel(XT),'KFold',10);
denoisedSegments = cell(size(XT));
for fold = 1:10
    fprintf('\n=============================\n');
    fprintf('Fold %d / 10\n',fold);
    fprintf('=============================\n');
    trainIdx = training(cv,fold);
    testIdx  = test(cv,fold);
    XTrain = XT(trainIdx);
    YTrain = YT(trainIdx);
    XTest = XT(testIdx);
    XTrain = cellfun(@single,XTrain,'UniformOutput',false);
    YTrain = cellfun(@single,YTrain,'UniformOutput',false);
    XTest = cellfun(@single,XTest,'UniformOutput',false);
    cvVal = cvpartition(numel(XTrain),'HoldOut',0.10);
    idxTrain = training(cvVal);
    idxVal   = test(cvVal);
    XVal = XTrain(idxVal);
    YVal = YTrain(idxVal);
    XTrain = XTrain(idxTrain);
    YTrain = YTrain(idxTrain);

    winLength = 32;
    overlapLength = 16;
    numFeatures=66;            
    win=rectwin(winLength);
    input = sequenceInputLayer(1, MinLength=640, Name="input");
    BN = batchNormalizationLayer(Name="BN");
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
    FC= fullyConnectedLayer(34, Name="fcMerge");
    FC2= fullyConnectedLayer(1, Name="fcFinal");
    Regression = regressionLayer(Name="regressionOutput");
    lgraph = layerGraph();
    lgraph = addLayers(lgraph, input);
    lgraph = addLayers(lgraph, BN);
    lgraph = addLayers(lgraph, [convA1; reluA1; convA2; reluA2]);          
    lgraph = addLayers(lgraph, [convB1; reluB1; transConvB2; reluTransB2]);
    lgraph = addLayers(lgraph, concatSkip);                                
    lgraph = addLayers(lgraph, [transConvC1; reluTransC1; transConvC2; reluTransC2; flattenC2; concat; bilstm; drop; FC;FC2;Regression]); 
    lgraph = connectLayers(lgraph, "input", "BN");
    lgraph = connectLayers(lgraph, "BN", "convA1");
    lgraph = connectLayers(lgraph, "reluA2", "convB1");
    lgraph = connectLayers(lgraph, "reluA2", "concatSkip/in1");
    lgraph = connectLayers(lgraph, "reluTransB2", "concatSkip/in2");
    lgraph = connectLayers(lgraph, "concatSkip", "transConvC1");
    flatten = flattenLayer(Name='FT');
    norm0 = layerNormalizationLayer(Name="norm0");
    attention = selfAttentionLayer(8,128,'AttentionMask','causal', Name = "ATT");
    fc1_1 = fullyConnectedLayer(256, Name="fc1_1");
    relu1 = geluLayer(Name="relu1");
    fc1_2 = fullyConnectedLayer(34, Name="fc1_2");
    add1 = additionLayer(2, Name="add1");
    norm1 = layerNormalizationLayer(Name="norm1");
    attention2 = selfAttentionLayer(16,64, Name='attention2');
    fc2_1 = fullyConnectedLayer(256, Name="fc2_1");
    relu2 = geluLayer(Name="relu2");
    fc2_2 = fullyConnectedLayer(34, Name="fc2_2");
    add2 = additionLayer(2, Name="add2");
    norm2 = layerNormalizationLayer(Name="norm2");
    attention3 = selfAttentionLayer(8,128, Name='attention3');
    fc3_1 = fullyConnectedLayer(256, Name="fc3_1");
    relu3 = geluLayer(Name="relu3");
    fc3_2 = fullyConnectedLayer(34, Name="fc3_2");
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
    lgraph = connectLayers(lgraph, "BN", "FT");
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
layers = lgraph;
analyzeNetwork(layers)
options = trainingOptions("adam", ...
    MaxEpochs=200, ...
    MiniBatchSize=16, ...
    InitialLearnRate=0.001, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.5, ...
    LearnRateDropPeriod=20, ...
    ValidationData={XVal,YVal},...
    ValidationFrequency=20, ...
    ValidationPatience=60, ...
    L2Regularization=0.001, ...
    GradientThreshold=1, ...
    Plots="training-progress", ...
    Shuffle="every-epoch", ...
    Verbose=false);
 net = trainNetwork(XTrain, YTrain, layers, options)
    % Predict
 YPred = predict(net, XTest, MiniBatchSize=16);
 denoisedSegments(testIdx) = YPred;
end

%  Reconstruction 
numSignals = length(detectedContSignals);
reconstructedOutput = cell(numSignals,1);
for i = 1:numSignals
    idx = find(segmentInfo(:,1) == i);
    [~, order] = sort(segmentInfo(idx,4));
    idx = idx(order);
    recSig = [];
    for k = 1:length(idx)
        seg = denoisedSegments{idx(k)};
        seg = seg(:);
        recSig = [recSig; seg];
    end
    reconstructedOutput{i} = recSig';
end
plot(reconstructedOutput{250, 1})
hold on
plot(detectedContSignals{250, 1})
hold on
plot(detectedCleanSignals{250, 1})
legend('pred','cont', 'clean')

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
plot(cont_theta_delta{30, 1}(22,:))
hold on
plot(reconstructedThetaDelta{30, 1}(22,:))
hold on
plot(clean_theta_delta{30, 1}(22,:))
legend('contaminated', 'predicted','clean')

cont_reconstructed = cont;   
for subj = 1:numSubjects
    for ch = 1:numChannels
        td = reconstructedThetaDelta{subj}(ch,:).';
        cont_reconstructed{subj}{ch}{4} = single(td);
        cont_reconstructed{subj}{ch}{5} = zeros(size(td), 'single');
    end
end
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
numSubjects = numel(clean);
numChannels = numel(clean{1});
clean_ = cell(numSubjects,1);
for subj = 1:numSubjects
    N = length(clean{subj}{1}{1});
    temp = zeros(numChannels, N, 'single');
    for ch = 1:numChannels
        gamma = clean{subj}{ch}{1};
        beta  = clean{subj}{ch}{2};
        alpha = clean{subj}{ch}{3};
        theta = clean{subj}{ch}{4};
        delta = clean{subj}{ch}{5};
        temp(ch,:) = (gamma + beta + alpha + theta + delta).';
    end
    clean_{subj} = temp;
end
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

plot(reconstructed_{20, 1}(23,1300:2000))
hold on
plot(clean_{20, 1}(23,1300:2000))
hold on
plot(cont_{20, 1}(23,1300:2000))
legend('pred', 'clean', 'cont')
