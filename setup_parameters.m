% Script to setup parameters used for the whole project

% parameters for the auto_win_select.m
parameters.largest_epidist_range = 3000;
parameters.cycle_before = 2;
parameters.cycle_after = 5;
parameters.min_dist_tol = deg2km(20);
parameters.max_dist_tol = deg2km(160);

% parameters for the gsdfmain.m
%
parameters.is_rm_resp = 1;
parameters.minstadist = 5;
parameters.maxstadist = 200;
parameters.periods = [20 25 32 40 50 60 80 100];
parameters.periods = sort(parameters.periods);  % make sure periods are ascending
parameters.refv = 4;
parameters.refphv = ones(size(parameters.periods))*4;
parameters.min_width = 0.06;
parameters.max_width = 0.10;
parameters.wintaperlength = 30;
parameters.prefilter = [10,200];
parameters.xcor_win_halflength = 100;
parameters.Nfit = 2;
parameters.Ncircle = 2;
parameters.cohere_tol = 0.5;
parameters.tp_tol = 10;

% parameters for the tomography
%
parameters.lalim=[-11.2 -7.8];
parameters.lolim=[148.8 151.5];
parameters.gridsize=0.1;
parameters.raydensetol=deg2km(parameters.gridsize)*2;
parameters.sou_dist_tol = 1;  % count by wavelength
parameters.smweight0 = 2.0;
parameters.maxerrweight = 2;
parameters.Tdumpweight = 0;
parameters.Rdumpweight = 0;
parameters.fiterrtol = 3;
parameters.dterrtol = 1;
parameters.isRsmooth = 0;
parameters.inverse_err_tol = 2;  % count be number of standard devition
parameters.min_amp_tol = 0.4;  % station with amplitude smaller than this ratio of average amplitude will not be used.


% parameter for stacking 
parameters.mincsnum=50;
parameters.min_phv_tol = 3;
parameters.max_phv_tol = 5;
parameters.is_raydense_weight = 0;
parameters.min_event_num = 10;
parameters.err_std_tol = 2;
