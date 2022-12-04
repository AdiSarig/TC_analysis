function [params] = init_params(load_masking, load_IB, load_DCF, debug)

params.debug = debug;

paradigm_names = {'Masking', 'Inattentional_Blindness', 'color_fusion'};
params.paradigm_names = paradigm_names([load_masking, load_IB, load_DCF]);

%% load data

% Set condition names and trigger values by condition
params.cond_names.ib = {'P1_f', 'P1_h', 'P1_n', 'P2_f', 'P2_h', 'P2_n', 'P3_f', 'P3_h', 'P3_n'};
params.cond_code.ib.P1_f   = 'S210'; params.cond_code.ib.P1_h   = 'S211'; params.cond_code.ib.P1_n   = 'S212';
params.cond_code.ib.P2_f   = 'S220'; params.cond_code.ib.P2_h   = 'S221'; params.cond_code.ib.P2_n   = 'S222';
params.cond_code.ib.P3_f   = 'S230'; params.cond_code.ib.P3_h   = 'S231'; params.cond_code.ib.P3_n   = 'S232';

params.cond_names.masking = {'V_f', 'V_h', 'V_b', 'IV_f', 'IV_h', 'IV_b'};
params.cond_code.masking.V_f  = 'S 10'; params.cond_code.masking.V_h  = 'S 12'; params.cond_code.masking.V_b  = 'S 14';
params.cond_code.masking.IV_f = 'S 11'; params.cond_code.masking.IV_h = 'S 13'; params.cond_code.masking.IV_b = 'S 15';

params.cond_names.dcf = {'V_f', 'V_h', 'V_b', 'IV_f', 'IV_h', 'IV_b'};
params.cond_code.dcf.V_f  = {'S 11', 'S 12'};
params.cond_code.dcf.V_h  = {'S 21', 'S 22'};
params.cond_code.dcf.V_b  = {'S 31', 'S 32'};
params.cond_code.dcf.IV_f = {'S 61', 'S 62'};
params.cond_code.dcf.IV_h = {'S 71', 'S 72'};
params.cond_code.dcf.IV_b = {'S 81', 'S 82'};

% load biosemi channel neighbours
load('biosemi64_neighb.mat');
params.neighbours = neighbours;

% define trial length
params.prestim = 1;
params.poststim = 1.5;

%% conditions, subjects, latencies & channels

% condition names without control
params.cond_names.masking_no_blank = {'V_f', 'V_h', 'IV_f', 'IV_h'};
params.cond_names.dcf_no_blank     = {'V_f', 'V_h', 'IV_f', 'IV_h'};
params.cond_names.ib_no_noise      = {'P1_f', 'P1_h', 'P2_f', 'P2_h'};
params.cond_names.ib_noise         = {'P1_n', 'P2_n'};

% select subjects
params.exclude_subjects.masking = [2, 4, 13]; % based on exclusion criteria
params.exclude_subjects.dcf     = [2, 3, 4, 13];
params.exclude_subjects.ib = [2, 3, 4, 13, 14];
params.exclude_subjects.ib_noticers = [1, 3, 5, 8, 10, 15, 1013];
% IBs = [2, 4, 6, 7, 9, 11, 12, 13, 14, 16, 17, 18, 1015, 1016];
params.TC_subs = [6:7, 9, 11:12, 16:18];

% select time points
params.p3_latency = [0.4, 0.6]; % based on Shafto & Pitts, 2015
params.van_latency.masking = [0.17, 0.25]; % based on visual inspection of the data
params.van_latency.dcf = [0.18, 0.32]; % based on visual inspection of the data
params.van_latency.ib = [0.14, 0.24]; % based on visual inspection of the data

params.van_latency.masking = [0.16, 0.19]; % based on reed lab
params.van_latency.dcf = [0.24, 0.32]; % based on reed lab
params.van_latency.ib = [0.17, 0.25]; % based on reed lab

% select channels
params.p3_channels = {'C3', 'C1', 'Cz', 'C2', 'C4',...
               'CP3', 'CP1', 'CPz', 'CP2', 'CP4',...
               'P3', 'P1', 'Pz', 'P2', 'P4'};
params.van_channels = {'P5', 'P7', 'PO3', 'PO7',...
                'O1', 'POz', 'Oz', 'Iz', 'O2'...
                'P6', 'P8', 'PO4', 'PO8'};
params.p3_channel_vec = zeros(length(params.p3_channels),1);
params.van_channel_vec = zeros(length(params.van_channels),1);

end

