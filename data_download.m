clear;

javaaddpath('IRIS-WS-2.0.6.jar');

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
component = parameters.component;

datacache = 'datacache';
if ~exist('eventmat','dir')
	mkdir('eventmat');
end

if exist('timestamp.mat','file') && parameters.is_use_timestamp
	temp = load('timestamp.mat');
	start_time = temp.start_time;
else
	start_time = parameters.start_time;
end
if exist('timestamp.mat','file') && ~parameters.is_use_timestamp
	delete('timestamp.mat');
end

network = parameters.network;
minMw = parameters.minMw;
maxdepth = parameters.maxdepth;

% end time is the current date minus two days to ensure the dateset completion
if isempty(parameters.end_time)
	end_time = datenum(date)-4;
	end_time = datestr(end_time,'yyyy-mm-dd HH:MM:SS');
else
	end_time = parameters.end_time;
end

events_info = irisFetch.Events('startTime',start_time,'endTime',end_time,...
		'MinimumMagnitude',minMw,'maximumDepth',maxdepth);

% data download and preparation
for ie=1:length(events_info)
	% gather event information
	otime = datenum(events_info(ie).PreferredTime,'yyyy-mm-dd HH:MM:SS.FFF');
	eventid = datestr(otime,'yyyymmddHHMM');
	event_filename = fullfile('eventmat',[eventid,'_',parameters.component,'.mat']);
	if exist(event_filename,'file')
		disp(['Exist: ',event_filename,', Skip!']);
		continue;
	end
	evla = events_info(ie).PreferredLatitude;
	evlo = events_info(ie).PreferredLongitude;
	evdp = events_info(ie).PreferredDepth;
	% init event structure and fill in event information 
	event.evla = evla;
	event.evlo = evlo;
	event.otime = otime*24*3600;
	event.id = eventid;
	event.otimestr = datestr(otime);

	% require station information
	event_starttime = datestr(otime,'yyyy-mm-dd HH:MM:SS');
	event_endtime = datestr(otime+parameters.datalength/3600/24,'yyyy-mm-dd HH:MM:SS');
	stations_info = irisFetch.Stations('channel',parameters.network,...
						'*','*',component,...
			'MinimumLatitude', lalim(1), 'MaximumLatitude', lalim(2), ...
			'MinimumLongitude',lolim(1) ,'MaximumLongitude',lolim(2), ...
						'startTime',event_starttime,'endTime',event_endtime);
	if ~exist(datacache,'dir')
		mkdir(datacache);
	end
	if ~exist(fullfile(datacache,eventid),'dir')
		mkdir(fullfile(datacache,eventid));
	end
	for ista =1:length(stations_info)
		stnm = stations_info(ista).StationCode;
		network = stations_info(ista).NetworkCode;
		sta_filename = fullfile('datacache',eventid,[eventid,'_',network,'_',stnm,'.mat']);
		if exist(sta_filename,'file')
			disp(['Exist: ',sta_filename,', Skip!']);
			continue;
		end
		disp(['Downloading station: ',stnm,' From:',event_starttime,' To:',event_endtime]);
		try
			if component(1) == 'L'
			traces = irisFetch.Traces(network,stnm,'*','LH*',event_starttime,event_endtime,'includePZ');
			elseif component(1) == 'B'
			traces = irisFetch.Traces(network,stnm,'*','BH*',event_starttime,event_endtime,'includePZ');
			end
			save(sta_filename,'traces');
		catch e
			e.message;
			continue;
		end
	end

	% remove the station lacking components
	sta_mat_files = dir(fullfile(datacache,eventid,'*.mat'));
	for ista = 1:length(sta_mat_files) 
		sta = load(fullfile(datacache,eventid,sta_mat_files(ista).name));
		ind = find(ismember({sta.traces.channel},{'LHZ','BHZ'}));
		if isempty(ind) || length(ind)>1 || sta.traces(ind).sampleCount < 500
			disp(['deleting ',sta_mat_files(ista).name]);
			delete(fullfile(datacache,eventid,sta_mat_files(ista).name));
			continue;
		end
		ind = find(ismember({sta.traces.channel},{'LHN','LH1','BHN','BH1'}));
		if isempty(ind) || length(ind)>1 || sta.traces(ind).sampleCount < 500
			disp(['deleting ',sta_mat_files(ista).name]);
			delete(fullfile(datacache,eventid,sta_mat_files(ista).name));
			continue;
		end
		ind = find(ismember({sta.traces.channel},{'LHE','LH2','BHE','BH2'}));
		if isempty(ind) || length(ind)>1 || sta.traces(ind).sampleCount < 500
			disp(['deleting ',sta_mat_files(ista).name]);
			delete(fullfile(datacache,eventid,sta_mat_files(ista).name));
			continue;
		end
	end
	
	% prepare the data, remove instrument response, transfer the format
	sta_mat_files = dir(fullfile(datacache,eventid,'*.mat'));
	if length(sta_mat_files) < 3
		continue;
	end
	for ista = 1:length(sta_mat_files)
		sta = load(fullfile(datacache,eventid,sta_mat_files(ista).name));
		ind = find(ismember({sta.traces.channel},{'LHZ','BHZ'}));
		LHZ = sta.traces(ind);
		ind = find(ismember({sta.traces.channel},{'LHN','LH1','BHN','BH1'}));
		LHN = sta.traces(ind);
		ind = find(ismember({sta.traces.channel},{'LHE','LH2','BHE','BH2'}));
		LHE = sta.traces(ind);
		stla = LHZ.latitude;
		stlo = LHZ.longitude;
		stnm = LHZ.station;
		stel = LHZ.elevation/1e3;
		if parameters.component == 'LHZ'
			LHZ = rm_resp(LHZ); 
			stadata = LHZ.data_cor;
			datacmp = 'LHZ';
			data_starttime = LHZ.startTime;
			data_delta = 1./LHZ.sampleRate;
		elseif parameters.component == 'BHZ'
			LHZ = rm_resp(LHZ); 
			stadata = LHZ.data_cor;
			datacmp = 'BHZ';
			data_starttime = LHZ.startTime;
			data_delta = 1./LHZ.sampleRate;
		elseif parameters.component == 'LHT'
			% do the rotation
			LHE = rm_resp(LHE); 
			LHN = rm_resp(LHN); 
			dataE = LHE.data_cor;
			dataN = LHN.data_cor;
			if LHE.sampleCount ~= LHZ.sampleCount
				% resample to make sure the two component have same number of data point
				b = (LHE.startTime - otime)*24*3600;
				delta = 1./LHE.sampleRate;
				E_timeaxis = b:delta:b + delta*(LHE.sampleCount-1);
				b = (LHN.startTime - otime)*24*3600;
				delta = 1./LHN.sampleRate;
				N_timeaxis = old_timeaxis(1):delta:old_timeaxis(end);
				dataN = interp1(N_timeaxis,dataN,E_timeaxis,'spline','extrap');
			end
			baz = azimuth(stla,stlo,evla,evlo);
			R_az = baz+180;
			T_az = R_az+90;
			dataR = dataN*cosd(R_az-LHN.azimuth)+dataE*cosd(R_az-LHE.azimuth);
			dataT = dataN*cosd(T_az-LHN.azimuth)+dataE*cosd(T_az-LHE.azimuth);
			stadata = dataT;
			datacmp = 'LHT';
			data_starttime = LHE.startTime;
			data_delta = 1./LHE.sampleRate;
		else
			disp('Unrecognized component, exit!');
			return;
		end % component selection
		% resample the data if necessary
		if parameters.resample_delta > data_delta
			old_taxis = 0:data_delta:(length(stadata)-1)*data_delta;
			resample_delta = parameters.resample_delta;
			new_taxis = 0:resample_delta:(length(stadata)-1)*data_delta;
			% apply anti-alias filter
			fN = 1/2/data_delta;
			w_c = 1./resample_delta/fN;
			[b,a] = butter(2,w_c,'low');
			stadata = filtfilt(b,a,stadata);
			stadata = interp1(old_taxis,stadata,new_taxis,'spline');
			data_delta = resample_delta;
		end
		event.stadata(ista).stla = stla;	
		event.stadata(ista).stlo = stlo;	
		event.stadata(ista).stel = stel;	
		event.stadata(ista).dist = vdist(stla,stlo,evla,evlo)/1e3;
		event.stadata(ista).otime = data_starttime*24*3600;
		event.stadata(ista).delta = data_delta;
		event.stadata(ista).data = stadata;	
		event.stadata(ista).cmp = datacmp;	
		event.stadata(ista).stnm = stnm;
	end  %ista
	save(event_filename,'event')
	disp(['Save to ',event_filename]);
end

if parameters.is_use_timestamp
	start_time = end_time;
	save('timestamp.mat','start_time');
end
