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
end
