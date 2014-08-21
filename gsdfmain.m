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

isdebug = 1;
is_overwrite = 0;

eventmatpath = './eventmat/';
CSoutputpath = './CSmeasure/';
winparapath = './winpara/';

if ~exist(CSoutputpath,'dir')
	mkdir(CSoutputpath)
end

% Setup parameters
setup_parameters

% Setup Error Codes for Bad data
setup_ErrorCode

comp = parameters.component;
periods = parameters.periods;
minstadist = parameters.minstadist;
maxstadist = parameters.maxstadist;
is_rm_resp = parameters.is_rm_resp;

respfile = 'staresp.mat';
if is_rm_resp
	temp = load(respfile);
	staresp = temp.staresp;
end

matfiles = dir([eventmatpath,'/*_',comp,'.mat']);
for ie = 1:length(matfiles)
%for ie = 100
    
	clear event eventcs CS
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;

	matfilename = [CSoutputpath,char(event.id),'_cs_',comp,'.mat'];
	if ~is_overwrite && exist(matfilename,'file')
		disp(['Found ',matfilename,', skip this event!']);
		continue;
	end
	if isfield(event,'isgood') && ~event.isgood
		disp(['Manually move',matfilename,', skip this event!']);
		continue;
	end
	disp(['Start to work on event: ',event.id]);

	% set up some useful arrays
	stlas = [event.stadata(:).stla];
	stlos = [event.stadata(:).stlo];
	stnms = {event.stadata(:).stnm};
	dists = [event.stadata(:).dist];

	% remove instrument response if necessary
	if is_rm_resp
		disp('Removing Station intrument response');
		resp_stanms = {staresp(:).staname};
		for ista = 1:length(event.stadata)
			resp_staid = find(ismember(resp_stanms,stnms(3)));
			respN = length(staresp(resp_staid).resp);
			if mod(respN,2) == 0
				resp_faxis = [0:respN/2,-respN/2+1:-1]/respN/staresp(resp_staid).dtr;  % only works for even data points
			else
				resp_faxis = [0:floor(respN/2),-floor(respN/2)+1:-1]/respN/staresp(resp_staid).dtr;  % only works for even data points
			end
			dataN = length(event.stadata(ista).data);
			if mod(dataN,2) == 0
				data_faxis = [0:dataN/2,-dataN/2+1:-1]/dataN/event.stadata(ista).delta;  
			else
				data_faxis = [0:floor(dataN/2),-floor(dataN/2):-1]/dataN/event.stadata(ista).delta;  
			end
			resp = interp1(resp_faxis,staresp(resp_staid).resp,data_faxis);
			fftdata = fft(event.stadata(ista).data);
			fftdata = fftdata.*resp(:);    % the resp has already been inversed in the input data.
%			event.stadata(ista).odata = event.stadata(ista).data;
			event.stadata(ista).data = real(ifft(fftdata));
		end
	end
	
	% find the win para files
	if ~isfield(event,'winpara') || length(event.winpara) == 1
		disp(['No win para information for event:',event.id]);
		continue;
	end

	% remove out of range stations
	for ista = 1:length(event.stadata)
		if event.stadata(ista).isgood > 0 &&...
			~Is_inrange(stlas(ista),stlos(ista),parameters)
			event.stadata(ista).isgood = ErrorCode.sta_outofrange;
		end
	end

	% Calculate Auto-correlation
	disp(['Calculating the auto-correlation of each station'])
	for ista = 1:length(event.stadata)
		if event.stadata(ista).isgood > 0
			autocor = CS_measure(event,ista,ista,parameters);
			if sum(autocor.amp) == 0
				disp(['Station: ',event.stadata(ista).stnm,' doesn''t have enough data for this event!']);
				event.stadata(ista).isgood = ErrorCode.sta_lackdata;
				event.autocor(ista) = autocor;
			else
				event.autocor(ista) = autocor;
			end
		else
			event.autocor(ista) = init_CSstruct;
		end
	end
	

	% start to find nearby stations and apply cross-correlation measurement
	csnum = 0;
	disp(['Calculating cross-correlation between stations'])
	for ista = 1:length(event.stadata)
		if event.stadata(ista).isgood < 0;
			continue;
		end
		% Find nearby stations
		stadist = deg2km(distance(stlas(ista),stlos(ista),stlas,stlos));
		nbstaids = find(stadist > minstadist & stadist < maxstadist);
		for nbsta = nbstaids
			if nbsta > ista && event.stadata(nbsta).isgood > 0
				% Build up Cross-Station Measurement structure
				csnum = csnum+1;
				if mod(csnum,100) == 0
					disp(csnum);
				end
				CS(csnum) = CS_measure(event,ista,nbsta,parameters);
			end % end of nbsta > ista
		end % end of nearby station loop
	end % end of station loop

	if csnum < 10
		disp(['Event has too few cross-correlations!']);
		continue;
	end
	% Calculate the coherency between each station pairs
	for ics = 1:length(CS)
		sta1 = CS(ics).sta1;
		sta2 = CS(ics).sta2;
		for ip=1:length(periods)
			CS(ics).cohere(ip) = CS(ics).amp(ip)^2/event.autocor(sta1).amp(ip)/event.autocor(sta2).amp(ip);
		end
	end

	% Removing window $W_S$ bias
	for ics = 1:length(CS)
		sta1 = CS(ics).sta1;
		sta2 = CS(ics).sta2;
		for ip=1:length(periods)
			CS(ics).dtp(ip) = CS(ics).dtp(ip)-event.autocor(sta2).dtp(ip);
		end
	end

	%% Data Quality Control
	%
	for ics = 1:length(CS)
		for ip = 1:length(periods)
			CS(ics).isgood(ip) = 1;
		end
	end

	% Remove measurements having fitting error
	for ics = 1:length(CS)
		for ip = 1:length(periods)
			if CS(ics).exitflag(ip) < 0
				CS(ics).isgood(ip) = ErrorCode.cs_fit_error;
			end
			if event.autocor(CS(ics).sta1).exitflag(ip) < 0
				CS(ics).isgood(ip) = ErrorCode.sta_fit_error;
			end
			if event.autocor(CS(ics).sta2).exitflag(ip) < 0
				CS(ics).isgood(ip) = ErrorCode.sta_fit_error;
			end
		end
	end
	% Remove measurements having low coherency
	for ics = 1:length(CS)
		for ip = 1:length(periods)
			if CS(ics).cohere(ip) < parameters.cohere_tol
				CS(ics).isgood(ip) = ErrorCode.low_cohere;
			end
		end
	end
	% Get average phase velocity across the array and remove the outliers.
	clear avgphv
	for ip=1:length(periods)
		clear ddist dtp isgood
		for ics = 1:length(CS)
			ddist(ics) = CS(ics).ddist;
			dtp(ics) = CS(ics).dtp(ip);
			isgood(ics) = CS(ics).isgood(ip);
		end % end of ics
		goodind = find(isgood > 0);
		para = polyfit(ddist(goodind),dtp(goodind),1);
		err = abs(ddist*para(1) + para(2) - dtp);
		for ics = 1:length(CS)
			if err(ics) > parameters.tp_tol && CS(ics).isgood(ip) > 0
				CS(ics).isgood(ip) = ErrorCode.high_tp_err;
			end
			isgood(ics) = CS(ics).isgood(ip);
		end
		goodind = find(isgood > 0);
		para = polyfit(ddist(goodind),dtp(goodind),1);
		avgphv(ip) = 1./para(1);
	end % end of periods

	% create eventcs structure and output
	eventcs.CS = CS;
	eventcs.autocor = event.autocor;
	eventcs.id = event.id;
	eventcs.avgphv = avgphv;
	eventcs.stlas = stlas;
	eventcs.stlos = stlos;
	eventcs.stnms = stnms;
	eventcs.evla = event.evla;
	eventcs.evlo = event.evlo;
	eventcs.dists = dists;
	eventcs.eventmatfile = [eventmatpath,matfiles(ie).name];

	matfilename = [CSoutputpath,char(event.id),'_cs_',comp,'.mat'];
	save(matfilename,'eventcs')
	disp(['Save to ',matfilename]);
end % end of ie loop
