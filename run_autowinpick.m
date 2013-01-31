% Scripts to run the auto_win_pick function for all the events and generate a old version "events" file
clear;

isdebug = 1;
is_overwrite = 0;

eventmatpath = './eventmat/';
outwinpath = './winpara/';

if ~exist(outwinpath,'dir')
	mkdir(outwinpath)
end

% Setup parameters
setup_parameters

% Setup Error Codes for Bad data
setup_ErrorCode

periods = parameters.periods;

matfiles = dir([eventmatpath,'/*.mat']);
for ie = 1:length(matfiles)

	clear event eventcs CS
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;
	disp(event.id)

	if ~is_overwrite
		filename = [outwinpath,'/',event.id,'.bad'];
		if exist(filename,'file')
			disp(['Exist: ',filename,' Skip!']);
			continue;
		end
		filename = [outwinpath,'/',event.id,'.win'];
		if exist(filename,'file')
			disp(['Exist: ',filename,' Skip!']);
			continue;
		end
	end

	% set up some useful arrays
	stlas = [event.stadata(:).stla];
	stlos = [event.stadata(:).stlo];
	stnms = {event.stadata(:).stnm};
	dists = [event.stadata(:).dist];

	% automatically select the signal window by using ftan method
	disp('Start to picking the window');
	tic
		[winpara event] = auto_win_select(event);
	toc
	if isdebug
		plot_win_select(event,periods,winpara);
	end
	if length(winpara) ~= 4
		filename = [outwinpath,'/',event.id,'.bad'];
		fp = fopen(filename,'w');
		fprintf(fp,'%f\n',0);
		fclose(fp);
		continue;
	end
	filename = [outwinpath,'/',event.id,'.win'];
	fp = fopen(filename,'w');
	fprintf(fp,'%s %f %f %f %f\n',event.id,winpara(1),winpara(2),winpara(3),winpara(4));
	fclose(fp);
	event.winpara = winpara;
	save([eventmatpath,matfiles(ie).name],'event');
	
end % end of event
