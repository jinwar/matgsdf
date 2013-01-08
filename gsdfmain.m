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

eventmatpath = './eventmat/';

minstadist = 5;
maxstadist = 200;
refv = 4;
periods = [20 25 32 40 50 60 80 100];


periods = sort(periods);  % make sure periods are ascending

matfiles = dir([eventmatpath,'/*.mat']);
%for ie = 1:length(matfiles)
for ie = 2:2
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

	% start doing the cross-correlation
	for ista = 1:length(event.stadata)
		stadist = deg2km(distance(stlas(ista),stlos(ista),stlas,stlos));
		nbstaids = find(stadist > minstadist & stadist < maxstadist);
		for nbsta = nbstaids
			if nbsta > ista


			end % end of nbsta > ista
		end % end of nearby station loop
	end % end of station loop
end % end of ie loop
