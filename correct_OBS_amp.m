%% script to correct the amplitude gain for OBSs. Should be run only once!!
% written by Ge Jin, jinwar@gmail.com
% March, 2013


OBS_cor = 3.0707
setup_parameters;
comp = parameters.component;


CSfiles = dir(['CSmeasure/*_',comp,'*.mat']);

OBSstnms = textread('OBS.lst','%s');

% Gather information
for ie = 1:length(CSfiles)
	clear eventcs amps
	load(fullfile('CSmeasure',CSfiles(ie).name));
	disp(CSfiles(ie).name)
	for ista = 1:length(eventcs.stnms)
		if ismember(eventcs.stnms(ista),OBSstnms)
%			disp(char(eventcs.stnms(ista)));
			eventcs.autocor(ista).amp = eventcs.autocor(ista).amp.*(OBS_cor.^2);
		end
	end
	save(fullfile('CSmeasure',CSfiles(ie).name),'eventcs');
end

