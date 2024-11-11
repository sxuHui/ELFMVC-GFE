function [LHs, LHLHs] = Hs2LHs(Hs, eta, knn_size,m)
nView = length(Hs);

% Ss = cell(1, nView);
% Ls = cell(1, nView);
% HLs = cell(1, nView);
LHs = cell(nView, nView);
for iView = 1:nView
    [~, Xa] = litekmeans(Hs{iView}, m, 'Replicates', 1);   
    B = ConstructBP_pkn(Hs{iView}, Xa, 'nNeighbor', knn_size);
    idx = sum(B, 1) > 0;
    B = B(:, idx);
    P = B * diag(sum(B, 1).^(-.5));
    w0 = exp(-eta) * eta^0/factorial(0);
    PTP = P' * P + 1e-5 * eye(size(P, 2));
    PTP = (PTP + PTP')/2;
    Lm = eye(size(PTP, 1)) - PTP;
    HL2 = expm(-eta * Lm); % m^3
    HL2 = (HL2 + HL2')/2;
    core_HK = (w0 * eye(size(B,2)) - HL2) / PTP; % m^3
    if sum(sum(isnan(core_HK))) + sum(sum(isinf(core_HK))) > 0
        disp('NAN');
        idx_filter = [idx_filter; iView1]; %#ok
    end
    for iView2 = 1:nView
        tmp1 = P' * Hs{iView2}; % m * n * d
        tmp2 = core_HK * tmp1; % m * m * d
        LHs{iView, iView2} = w0 * Hs{iView2} - P * tmp2; % n m d
    end
end

clear S DS L HL;
if nargout > 1
    nDim = size(LHs{1,1}, 2);
    LHLHs = cell(nView, 1);
    for iView1 = 1:nView
        
        idx = 0;
        ABs = zeros((nDim*nDim), nView*(nView+1)/2);
        for iView2 = 1:nView
            LHa = LHs{iView1, iView2};
            for iView3 = iView2:nView
                LHb = LHs{iView1, iView3};
                idx = idx + 1;
                LHaLHb = LHa' * LHb;
                ABs(:, idx) = LHaLHb(:);
            end
        end
        LHLHs{iView1} = ABs;
    end
end

end