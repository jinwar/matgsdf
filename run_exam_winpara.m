clear;

eventmat_files = dir('eventmat/*.mat');

for ie=1:length(eventmat_files)
	load(fullfile('eventmat',eventmat_files(ie).name));
	if isfield(event,'winpara')
		if ~isfield(event,'isgood')
			event.isgood = 1;
		end
		disp(eventmat_files(ie).name);
		disp('p for pick window, g for isgood, q for next event');
		event = exam_winpara(event);
	end
	save(fullfile('eventmat',eventmat_files(ie).name),'event');
end

