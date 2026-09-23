function threshold = computeMFSThreshold1(MFS, method)

switch lower(method)

    %==============================================================
    % Unsupervised - Otsu
    %==============================================================
    case 'otsu'

    % Require some spread in the scores
    if std(MFS) < 0.05
        threshold = inf;
        return;
    end

    % First Otsu threshold
    MFSn = rescale(MFS);
    level = graythresh(MFSn);
    threshold = min(MFS) + level*(max(MFS)-min(MFS));

    pred = MFS >= threshold;

    % If Otsu isolates only one channel, remove that channel and run Otsu again
    if sum(pred) == 1

        idxOut = find(pred);              % extreme high-value channel
        MFS2 = MFS;
        MFS2(idxOut) = [];                % remove it temporarily

        MFS2n = rescale(MFS2);
        level2 = graythresh(MFS2n);
        threshold2 = min(MFS2) + level2*(max(MFS2)-min(MFS2));

        % Use the second threshold
        threshold = threshold2;
    end

    % Final check
    pred = MFS >= threshold;

    %==============================================================
    % Unsupervised - Gaussian Mixture Model
    %==============================================================
    %==============================================================
% Unsupervised - Gaussian Mixture Model
%==============================================================
case 'gmm'

    % Very small spread -> all channels assumed clean
    if std(MFS) < 0.05
        threshold = inf;
        return;
    end

    try

        rng(1)

        options = statset('MaxIter',500);

        gm = fitgmdist(MFS,2,...
            'Replicates',10,...
            'RegularizationValue',1e-3,...
            'CovarianceType','diagonal',...
            'Options',options);

        % Means
        mu = gm.mu(:);

        % Variances
        sigma = squeeze(gm.Sigma);

        % Mixing proportions
        p = gm.ComponentProportion(:);

        % Sort components
        [mu,idx] = sort(mu);
        sigma = sigma(idx);
        p = p(idx);

        % Separation
        separation = abs(mu(2)-mu(1))/sqrt(mean(sigma));

        % Nearly one cluster
        if separation < 1
            threshold = inf;
            return;
        end

        % Mid-point approximation
        threshold = mean(mu);

        % Reject if almost every channel falls into one class
        pred = MFS >= threshold;

        if sum(pred)<2 || sum(~pred)<2
            threshold = inf;
        end

    catch

        % GMM failed
        threshold = inf;

    end
end