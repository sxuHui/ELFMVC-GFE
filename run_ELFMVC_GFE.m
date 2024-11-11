%
%
%

clear;
clc;
data_path = fullfile(pwd,  filesep, "data_Hs", filesep);
addpath(data_path);
lib_path = fullfile(pwd, filesep, "lib", filesep);
addpath(lib_path);

dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};


exp_n = 'ELFMVC_GFE';
% profile off;
% profile on;
for i1 = 1:length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    dir_name = [pwd, filesep, exp_n, filesep, data_name, filesep];
    create_dir(dir_name);
    fname2 = fullfile(data_path, [data_name, '.mat']);
    load(fname2);
    nCluster = length(unique(Y));
    nView = length(Hs);
    nDim = nCluster * 4;
    
    % embeddings_s = [1, 2, 3, 4];
    embeddings_s = [1]; % default
    % eta_s = [0, 1, 3, 5, 7, 9]; % default
    eta_s = [9]; % default
    % eta_s = [9]; % default
    % knn_size_s = [5, 10]; % default
    knn_size_s = [5]; % default
    m_s = nCluster * [4:2:10];
    % lambda_s = [10^(-2), 10^(-1), 10^(0), 10^(1), 10^(2)];
    lambda_s = [0:0.1:1];
    % lambda_s = [0.8];
    paramCell = cell(1, length(embeddings_s) * length(eta_s) * length(knn_size_s) * length(m_s) * length(lambda_s));
    idx = 0;
    for iParam1 = 1:length(embeddings_s)
        for iParam2 = 1:length(eta_s)
            for iParam3 = 1:length(knn_size_s)
                for iParam4 = 1: length(m_s)
                    for iParam5 = 1:length(lambda_s)
                        idx = idx + 1;
                        param = [];
                        param.nEmbedding = embeddings_s(iParam1);
                        param.eta = eta_s(iParam2);
                        param.knn_size = knn_size_s(iParam3);
                        param.m = m_s(iParam4);
                        param.lambda = lambda_s(iParam5);
                        paramCell{idx} = param;
                    end
                end
            end
        end
    end
    paramCell = paramCell(~cellfun(@isempty, paramCell));
    nParam = length(paramCell);
    
    nMeasure = 13;
    nRepeat = 10;
    seed = 2024;
    rng(seed);
    % Generate 50 random seeds
    random_seeds = randi([0, 1000000], 1, nRepeat);
    % Store the original state of the random number generator
    original_rng_state = rng;
    %*********************************************************************
    % ELFMVC_GFE
    %*********************************************************************
    fname2 = fullfile(dir_name, [data_name, '_ELFMVC_GFE.mat']);
    if ~exist(fname2, 'file')
        ELFMVC_GFE_result = zeros(nParam, 1, nRepeat, nMeasure);
        ELFMVC_GFE_time = zeros(nParam, 1);
        for iParam = 1:nParam
            param = paramCell{iParam};
            nEmbedding = param.nEmbedding * nCluster;
            eta = param.eta;
            knn_size = param.knn_size;
            m = param.m;
            lambda = param.lambda;
            Hs_new = cell(1, nView);
            for iKernel = 1:nView
                Hi = Hs{iKernel};
                Hs_new{iKernel} = Hi(:, 1: nEmbedding);
            end
            t1_s = tic;
            if eta > 0
                LHs = Hs2LHs(Hs_new, eta, knn_size, m);
            else
                LHs = Hs_new;
            end
            t1 = toc(t1_s);
            t2_s = tic;
            for iRepeat = 1:nRepeat
                disp(['Param ', num2str(iParam), ' Repeat ', num2str(iRepeat)]);
                % Restore the original state of the random number generator
                rng(original_rng_state);
                % Set the seed for the current iteration
                rng(random_seeds(iRepeat));
                Ha = cell2mat(Hs_new);
                Ha = bsxfun(@rdivide, Ha, sqrt(sum(Ha.^2, 2)) + eps);
                label_0 = litekmeans(Ha, nCluster, 'MaxIter', 50, 'Replicates', 1);
                Y_0 = ind2vec(label_0')';
                [Iabel, Ws, alpha, beta, objHistory] = ELFMVC_GFE(Hs_new, nCluster, LHs, Y_0, lambda);
                result_aio = my_eval_y(Iabel, Y);
                ELFMVC_GFE_result(iParam, 1, iRepeat, :) = result_aio';
            end
            t2 = toc(t2_s);
            ELFMVC_GFE_time(iParam) = t1 + t2/nRepeat;
        end
        a1 = sum(ELFMVC_GFE_result, 2);
        a3 = sum(a1, 3);
        a4 = reshape(a3, nParam, nMeasure);
        a4 = a4/nRepeat;
        ELFMVC_GFE_result_grid = a4;
        ELFMVC_GFE_result_summary = [max(a4, [], 1), sum(ELFMVC_GFE_time)];
        save(fname2, 'ELFMVC_GFE_result', 'ELFMVC_GFE_time', 'ELFMVC_GFE_result_summary', 'ELFMVC_GFE_result_grid');
        
        
        disp([data_name, ' has been completed!']);
    end
end
% profile viewer;
rmpath(data_path);
rmpath(lib_path);