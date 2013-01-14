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
for ie = 1
    
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

	% start doing the cross-correlation
	for ista = 1:length(event.stadata)

		% Find nearby stations
		stadist = deg2km(distance(stlas(ista),stlos(ista),stlas,stlos));
		nbstaids = find(stadist > minstadist & stadist < maxstadist);
		for nbsta = nbstaids
			if nbsta > ista
				% read in data for station 1
				data1 = event.stadata(ista).data;
				bgtime = event.stadata(ista).otime - event.otime;
				dt1 = event.stadata(ista).delta;
				Nt = length(event.stadata(ista).data);
				taxis1 = bgtime + [0:Nt-1]'*dt1;
				% read in data for station 2
				data2 = event.stadata(nbsta).data;
				bgtime = event.stadata(nbsta).otime - event.otime;
				dt2 = event.stadata(nbsta).delta;
				Nt = length(event.stadata(nbsta).data);
				taxis2 = bgtime + [0:Nt-1]'*dt2;
				% resample the data if necessary
				if dt1 > dt2
					new_taxis2 = taxis2(1):dt1:taxis2(end);
					data2 = interp1(taxis2,data2,new_taxis2);
					taxis2 = new_taxis2;
				elseif dt1 < dt2
					new_taxis1 = taxis1(1):dt2:taxis1(end);
					data1 = interp1(taxis1,data1,new_taxis1);
					taxis1 = new_taxis1;
				end
			end % end of nbsta > ista
		end % end of nearby station loop
	end % end of station loop
end % end of ie loop
