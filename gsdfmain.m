% This script is used to measure the phase difference between nearby stations. 
% The input format is eventmat format, which is a matlab structure includes:
% struct event:
%   * string: eventid
%   * string: otimestr
%   * float otime
%   * float evla
%   * float evlo
%   * struct[] stadata
%      * stla
%      * stlo
%      * staotime
%      * delta
%      * data
%      * string comp

clear;

isdebug = 0;

eventmatpath = './eventmat/';

parameters.minstadist = 5;
parameters.maxstadist = 200;
parameters.refv = 4;
parameters.periods = [20 25 32 40 50 60 80 100];
parameters.min_width = 0.06;
parameters.max_width = 0.10;
parameters.wintaperlength = 30;
parameters.prefilter = [10,200];
parameters.xcor_win_halflength = 100;
parameters.Nfit = 2;

parameters.periods = sort(parameters.periods);  % make sure periods are ascending

periods = parameters.periods;
minstadist = parameters.minstadist;
maxstadist = parameters.maxstadist;

matfiles = dir([eventmatpath,'/*.mat']);
%for ie = 1:length(matfiles)
for ie = 100
    
	clear event csmeasure
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;

	% set up some useful arrays
	stlas = [event.stadata(:).stla];
	stlos = [event.stadata(:).stlo];
	dists = [event.stadata(:).dist];

	% automatically select the signal window by using ftan method
	winpara = auto_win_select(event,periods);
    plot_win_select(event,periods,winpara);
	if length(winpara) ~= 4
		continue;
	end
	v1 = winpara(1); t1=winpara(2); v2=winpara(3); t2=winpara(4);
	event.winpara = winpara;

	% start doing the cross-correlation
	csnum = 0;
	for ista = 1:length(event.stadata)

		% Find nearby stations
		stadist = deg2km(distance(stlas(ista),stlos(ista),stlas,stlos));
		nbstaids = find(stadist > minstadist & stadist < maxstadist);
		for nbsta = nbstaids
			if nbsta > ista
				% Build up Cross-Station Measurement structure
				csnum = csnum+1;
				disp(csnum);
				CS(csnum) = CS_measure(event,ista,nbsta,parameters);
			end % end of nbsta > ista
		end % end of nearby station loop
	end % end of station loop
end % end of ie loop
