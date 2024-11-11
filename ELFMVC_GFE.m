function [Iabel, Ws, alpha, beta, objHistory] = ELFMVC_GFE(Hs, nCluster, LHs, Y, lambda, Y_0)
nView = length(Hs);
nSmp = size(Hs{1}, 1);

%*************************************************
% Initialization of Ws
%*************************************************
Ws = cell(1, nView);
for iView = 1 : nView
    Ws{iView} = eye(size(Hs{iView}, 2), nCluster);
end

%*************************************************
% Initialization of alpha, beta
%*************************************************
alpha = ones(nView, 1)*nView;
beta = ones(nView, 1)/nView;
%*************************************************
% Initialization of Y
%*************************************************
if ~exist('Y', 'var')
    Ha = cell2mat(Hs);
    Ha = bsxfun(@rdivide, Ha, sqrt(sum(Ha.^2, 2)) + eps);
    label = litekmeans(Ha, nCluster, 'MaxIter', 50, 'Replicates', 10);
    Y = ind2vec(label')';
    clear Ha;
end

betaLHs = cell(1, nView);
for iView1 = 1:nView
    betaLH = zeros(nSmp, size(Hs{iView1}, 2));
    for iView2 = 1:nView
        betaLH = betaLH + beta(iView2) * LHs{iView1, iView2};
    end
    betaLHs{iView1} = betaLH;
end

iter = 0;
objHistory = [];
% accHistory = [];
converges = false;
maxIter = 30;
while ~converges
    %***********************************************
    % Update Ws
    %***********************************************
    for iView = 1 : nView
        HLLHS = (1 - lambda) * Hs{iView}' + lambda * betaLHs{iView}';
        HLLHSY = HLLHS * Y;
        [U, ~, V] = svd(HLLHSY, 'econ');
        Ws{iView} = U * V';
    end
    % [obj2, o1, o2] = compute_obj(Hs, LHs, Y, Ws, lambda, alpha, beta);
    % objHistory = [objHistory; obj2]; %#ok
    
    %***********************************************
    % Update Y
    %***********************************************
    HLLHW = zeros(nSmp, nCluster);
    for iView = 1:nView
        HLLHS = (1 - lambda) * Hs{iView} + lambda * betaLHs{iView};
        HLLHW = HLLHW + alpha(iView) * HLLHS * Ws{iView}; % n d c m
    end
    [~, label] = max(HLLHW, [], 2);
    Y = full(sparse(1:nSmp, label, ones(nSmp, 1), nSmp, nCluster));
    % [obj2, o1, o2] = compute_obj(Hs, LHs, Y, Ws, lambda, alpha, beta);
    % objHistory = [objHistory; obj2]; %#ok
    
    %***********************************************
    % Update beta
    %***********************************************
    A = zeros(nView, nView);
    b = zeros(nView, 1);
    for iView1 = 1 : nView
        Ai = zeros(nView, nView);
        bi = zeros(nView, 1);
        LHWs = cell(1, nView);
        for iView2 = 1 : nView
            LHWs{iView2} = LHs{iView1, iView2} * Ws{iView1}; % n d c m^2
        end
        for iView2 = 1 : nView
            for iView3 = iView2 : nView
                Ai(iView2, iView3) = sum(sum( LHs{iView1, iView2}.* LHs{iView1, iView3})); % n c m^3
            end
            LHWs = LHs{iView1, iView2} * Ws{iView1};
            bi(iView2) = sum(sum( LHWs .* Y));
        end
        Ai = max(Ai, Ai');
        A = A + alpha(iView1) * Ai;
        b = b + alpha(iView1) * bi;
    end
    opt = [];
    opt.Display = 'off';
    beta_old = beta;
    [beta, fval, ~, ~] = quadprog(A, -b, [], [], ones(1, nView), 1, zeros(nView, 1), ones(nView, 1), beta_old, opt);
    % obj = 2 * fval;
    % [obj2, o1, o2] = compute_obj(Hs, LHs, Y, Ws, lambda, alpha, beta);
    % objHistory = [objHistory; obj2]; %#ok

    betaLHs = cell(1, nView);
    for iView = 1:nView
        betaLH = zeros(nSmp, size(Hs{iView}, 2));
        for iView2 = 1:nView
            betaLH = betaLH + beta(iView2) * LHs{iView, iView2};
        end
        betaLHs{iView} = betaLH;
    end

    % obj2 = 0;
    % for iView = 1:nView
    %   o1 = sum(sum(betaLHs{iView} .* betaLHs{iView}));
    %   WH = betaLHs{iView} * Ws{iView};
    %   o2 = sum(sum(WH .* Y));
    %   obj2 = obj2 + alpha(iView) * (o1 - 2 * o2); 
    % end
    
    %***********************************************
    % Update alpha
    %***********************************************
    es = zeros(nView, 1);
    for iView = 1 : nView
        YWs = Y * Ws{iView}';
        previous = Hs{iView} - YWs; % n d c m
        latter = betaLHs{iView} - YWs;
        es(iView) = (1 - lambda) * sum(sum( previous.^2 )) + lambda * sum(sum( latter.^2 ));
    end
    alpha = sum(sqrt(es))./sqrt(es);
    obj = sum(alpha .* es);
    % obj = sum(es.*alpha);
    % [obj2, o1, o2] = compute_obj(Hs, LHs, Y, Ws, lambda, alpha, beta);

    objHistory = [objHistory; obj]; %#ok
    % [~, Iabel] = max(Y, [], 2);
    % result_aio = my_eval_y(Iabel, Y_0);
    % accHistory = [accHistory, result_aio(1)]; %#ok
    
    iter = iter + 1;
    if iter >= maxIter
        converges = 1;
    end
end
[~, Iabel] = max(Y, [], 2);
end

function [obj, o1, o2] = compute_obj(Hs, LHs, Y, Ws, lambda, alpha, beta)
nView = length(Hs);
nSmp = size(Hs{1}, 1);

betaLHs = cell(1, nView);
for iView1 = 1:nView
    betaLH = zeros(nSmp, size(Hs{iView1}, 2));
    for iView2 = 1:nView
        betaLH = betaLH + beta(iView2) * LHs{iView1, iView2};
    end
    betaLHs{iView1} = betaLH;
end

previous = zeros(nView, 1);
latter = zeros(nView, 1);
for iView = 1:nView
    pre = Hs{iView} - Y * Ws{iView}';
    previous(iView) = sum(sum(pre .* pre));
    lat = betaLHs{iView} - Y * Ws{iView}';
    latter(iView) = sum(sum(lat .* lat));
end
o1 = sum(alpha .* previous);
o2 = sum(alpha .* latter);
obj = (1 - lambda) * o1 + lambda * o2;
end