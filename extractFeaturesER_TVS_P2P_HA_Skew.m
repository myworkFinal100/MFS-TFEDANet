function Features = extractFeaturesER_TVS_P2P_HA_Skew(theta_delta, bands)
%=========================================================
% theta_delta : cell array
%               Each cell = Channels x Samples
%
% bands : FDM decomposition
%
% Output:
% Features : (Subjects*Channels) x 5
%
% Columns
% 1. ER_delta
% 2. TVS
% 3. Peak-to-Peak (P2P)
% 4. Hjorth Activity
% 5. Skewness
%=========================================================

numSubjects = numel(theta_delta);
numChannels = size(theta_delta{1},1);

Features = zeros(numSubjects*numChannels,5);

row = 1;

for subj = 1:numSubjects

    for ch = 1:numChannels

        %%=========================================================
        % Theta + Delta signal
        %%=========================================================
        sig = double(theta_delta{subj}(ch,:));

        %%=========================================================
        % Feature 1 : Delta Energy Ratio
        %%=========================================================
        delta = double(bands{subj}{ch}{5});
        theta = double(bands{subj}{ch}{4});
        alpha = double(bands{subj}{ch}{3});
        beta  = double(bands{subj}{ch}{2});
        gamma = double(bands{subj}{ch}{1});

        E_delta = sum(delta.^2);

        E_total = E_delta + ...
                  sum(theta.^2) + ...
                  sum(alpha.^2) + ...
                  sum(beta.^2)  + ...
                  sum(gamma.^2);

        ER_delta = E_delta/(E_total + eps);

        %%=========================================================
        % Feature 2 : Temporal Variability Score (TVS)
        %%=========================================================
        dx = diff(sig);

        TVS = std(dx)/(median(abs(dx)) + eps);

        %%=========================================================
        % Feature 3 : Peak-to-Peak Amplitude (P2P)
        %%=========================================================
        P2P = max(sig) - min(sig);

        %%=========================================================
        % Feature 4 : Hjorth Activity
        %%=========================================================
        HjorthAct = var(sig);

        %%=========================================================
        % Feature 5 : Skewness
        %%=========================================================
        Skew = skewness(sig);

        %%=========================================================
        % Store Features
        %%=========================================================
        Features(row,:) = [ER_delta TVS P2P HjorthAct Skew];

        row = row + 1;

    end

end

end