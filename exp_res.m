clear;
clc;
data_path = fullfile(pwd,  filesep, "data_Hs", filesep);
addpath(data_path);
lib_path = fullfile(pwd,  filesep, "lib", filesep);
addpath(lib_path);

dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};
datasetCandi = datasetCandi';
exp_s = {'avgH',... 
    'AWP',... % 2018
    'LFMVC',... % 2019
    'OPLF',... % 2021
    'ALMVC',... % 2022
    'LFLKA',... % 2023
    'M3LF', ... % 2024
    'ERMKC', ... % 2024
    'sLGm',... % 2024
    'HKLFMVC_discrete_global',... % 2024
    'XXXX',...
    'avg_GH',...
    'LFKMMVC',...
    'ELFMVC_GFE',...
    };
% exp_s = {'HKLFMVCNF',...
%     };


grid_size = zeros(1, length(exp_s));
smp_size = zeros(1, length(datasetCandi));
res_aio_cell = cell(length(datasetCandi), length(exp_s));

acc_table = zeros(length(exp_s), length(datasetCandi));
nmi_table = zeros(length(exp_s), length(datasetCandi));
ari_table = zeros(length(exp_s), length(datasetCandi));
f1_table = zeros(length(exp_s), length(datasetCandi));
% pur_table = zeros(length(exp_s), length(datasetCandi));
% ri_table = zeros(length(exp_s), length(datasetCandi));

acc_std_table = zeros(length(exp_s), length(datasetCandi));
nmi_std_table = zeros(length(exp_s), length(datasetCandi));
ari_std_table = zeros(length(exp_s), length(datasetCandi));
f1_std_table = zeros(length(exp_s), length(datasetCandi));
% pur_std_table = zeros(length(exp_s), length(datasetCandi));
% ri_std_table = zeros(length(exp_s), length(datasetCandi));

time_table = zeros(length(exp_s), length(datasetCandi));

for i1 =1: length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    
    % sample_size_str = regexp(data_name, '\d+n', 'match'); % Extract the substring containing the sample size
    % smp_size(i1) = str2double(sample_size_str{1}(1:end-1));
    
    res_i = zeros(length(exp_s), 11);
    for i2 =1: length(exp_s)
        exp_n = exp_s{i2};
        dir_name = [pwd, filesep, exp_n, filesep, data_name];
        fname2 = fullfile(dir_name, [data_name, '_', exp_n, '.mat']);
        
        if exist(fname2, 'file')
            load(fname2);
            s1 = [exp_n, '_avg_result_grid = reshape(mean(mean(', exp_n, '_result,3), 2), size(', exp_n, '_result, 1), size(', exp_n, '_result, 4));'];
            eval(s1);
            % s2 = ['[~, ', exp_n, '_bst_param_idx] = max(sum(', exp_n, '_avg_result_grid(:, [1,2,4,8]), 2));'];
            % s2 = ['[~, ', exp_n, '_bst_param_idx] = max(sum(', exp_n, '_avg_result_grid(:, [1,2,3,5]), 2));'];
            s2 = ['[~, ', exp_n, '_bst_param_idx] = max(sum(', exp_n, '_avg_result_grid(:, [1]), 2));'];
            eval(s2);
            s3 = [exp_n, '_result_selected = ', exp_n, '_result(', exp_n, '_bst_param_idx, :, :, :);'];
            eval(s3);
            s4 = [exp_n, '_time_selected = ', exp_n, '_time(', exp_n, '_bst_param_idx);'];
            eval(s4);
            s5 = [exp_n, '_result_selected_avg = mean(reshape(', exp_n, '_result_selected, size(', exp_n, '_result_selected, 3) * size(', exp_n, '_result_selected, 2), size(', exp_n, '_result_selected, 4)));'];
            eval(s5);
            s6 = [exp_n, '_result_selected_std = std(reshape(', exp_n, '_result_selected, size(', exp_n, '_result_selected, 3) * size(', exp_n, '_result_selected, 2), size(', exp_n, '_result_selected, 4)));'];
            eval(s6);
            s7 = [exp_n, '_result_selected_avg_time = ', exp_n, '_time_selected;'];
            eval(s7);

            s8 = ['acc_table(i2, i1) = ', exp_n, '_result_selected_avg(1);'];
            eval(s8); %#ok
            s9 = ['nmi_table(i2, i1) = ', exp_n, '_result_selected_avg(2);'];
            eval(s9); %#ok
            s10 = ['ari_table(i2, i1) = ', exp_n, '_result_selected_avg(4);'];
            % s10 = ['pur_table(i2, i1) = ', exp_n, '_result_selected_avg(3);'];
            eval(s10); %#ok
            s11 = ['f1_table(i2, i1) = ', exp_n, '_result_selected_avg(8);'];
            % s11 = ['ri_table(i2, i1) = ', exp_n, '_result_selected_avg(5);'];
            eval(s11); %#ok

            s12 = ['acc_std_table(i2, i1) = ', exp_n, '_result_selected_std(1);'];
            eval(s12); %#ok
            s13 = ['nmi_std_table(i2, i1) = ', exp_n, '_result_selected_std(2);'];
            eval(s13); %#ok
            s14 = ['ari_std_table(i2, i1) = ', exp_n, '_result_selected_std(4);'];
            % s14 = ['pur_std_table(i2, i1) = ', exp_n, '_result_selected_std(3);'];
            eval(s14); %#ok
            s15 = ['f1_std_table(i2, i1) = ', exp_n, '_result_selected_std(8);'];
            % s15 = ['ri_std_table(i2, i1) = ', exp_n, '_result_selected_std(5);'];
            eval(s15); %#ok

            s16 = ['time_table(i2, i1) = ', exp_n, '_result_selected_avg_time(1);']; 
            eval(s16); %#ok

            s17 = ['clear ', exp_n, '_avg_result_grid ', exp_n, '_bst_param_idx ', exp_n, '_result_selected ', exp_n, '_time_selected'];
            eval(s17); 
            s18 = ['clear ', exp_n, '_result_selected_avg ', exp_n, '_result_selected_std ', exp_n, '_result_selected_avg_time'];
            eval(s18);
        end
    end
end

