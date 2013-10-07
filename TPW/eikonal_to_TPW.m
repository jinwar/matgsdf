%% This script is used for convert the output of matGSDF into input files for Two-Plane-Wave fortran
%  code. 
%  Written by Ge Jin, jinwar@gmail.com
%  Sep 2013

clear;

min_sta_num = 10;

setup_parameters;
comp = parameters.component;
periods = parameters.periods;
refvs = [3.53 3.61 3.7 3.76 3.82 3.85];


% input path
eventcs_path = './CSmeasure/';
eikonal_output_path = './eikonal/';

for ip = 1:6
% output files
period = periods(ip);
control_file = ['all.0',num2str(period)];
ampph_file = ['phampcor.0',num2str(period),'.inp'];

% input parameters
OBS_amp_bust = 3.07;
nfreq = 1;
freq = 1/periods(ip);
detailoutput = ['detail.0',num2str(period),'.inp'];
summaroutput = ['summar.0',num2str(period),'.inp'];
gridnodes = 'PNG_gridnodes';
varfile = ['covar.0',num2str(period),'.inp'];
mavamp = ['mavamp.0',num2str(period),'.inp'];
fvelarea0 = ['temp.0',num2str(period),'.inp'];
velarea = ['velarea.0',num2str(period),'.inp'];
iterlimit = 10;
wlambda = 040;
dampvel = 0.25;
dampaniso = 0.01;
unifvel = refvs(ip);
kernelfiledir = dir(['kernel.0',num2str(period),'s*']);
kernelfile = kernelfiledir(1).name;

% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	disp('Found Bad stations:')
	disp(badstnms)
end

% gather necessary information from GSDF output
eikonal_matfiles = dir([eikonal_output_path,'/*eikonal_',comp,'.mat']);
goodeventnum = 0;
for ie=1:length(eikonal_matfiles)
	temp = load([eikonal_output_path,eikonal_matfiles(ie).name]);
	eventphv = temp.eventphv;
	disp(eventphv(ip).id);
	stlas = eventphv.stlas;
	stlos = eventphv.stlos;
	stnms = eventphv.stnms;
	if exist('badstnms','var')
		badstaids = find(ismember(stnms,badstnms));
	else
		badstaids = [];
	end
	tp = eventphv(ip).traveltime;
	goodstaind = find(~isnan(tp));
	if ~isempty(badstaids)
		badind = find(ismember(goodstaind,badstaids));
		goodstaind(badind) = [];
	end
	goodstanum = length(goodstaind);
	if goodstanum > min_sta_num
		temp = load([eventcs_path,'/',eventphv(ip).id,'_cs_',comp,'.mat']);
		eventcs = temp.eventcs;
		goodeventnum = goodeventnum + 1;
		twpevents(goodeventnum).evla = eventphv(ip).evla;
		twpevents(goodeventnum).evlo = eventphv(ip).evlo;
		twpevents(goodeventnum).stlas = eventphv(ip).stlas(goodstaind);
		twpevents(goodeventnum).stlos = eventphv(ip).stlos(goodstaind);
		twpevents(goodeventnum).stnms = eventphv(ip).stnms(goodstaind);
		twpevents(goodeventnum).tp = eventphv(ip).traveltime(goodstaind);
		twpevents(goodeventnum).id = eventphv(ip).id;
		for i = 1:length(goodstaind)
			ista = goodstaind(i);
			twpevents(goodeventnum).amp(i) = eventcs.autocor(ista).amp(ip);
		end
	end
end

% Calculate useful informations for TPW output
for ie=1:length(twpevents)
	stlas = twpevents(ie).stlas;
	stlos = twpevents(ie).stlos;
	evla = twpevents(ie).evla;
	evlo = twpevents(ie).evlo;
	degs = distance(stlas,stlos,evla,evlo);
	twpevents(ie).degs = degs;
	twpevents(ie).dists = deg2km(degs);
	twpevents(ie).azi = azimuth(evla,evlo,stlas,stlos);
	twpevents(ie).baz = azimuth(stlas,stlos,evla,evlo);
end

% output the files
cfp = fopen(control_file,'w');
afp = fopen(ampph_file,'w');

fprintf(cfp,'%d\n',length(twpevents));
for ie=1:length(twpevents)
	fprintf(cfp,'%d   %d\n',length(twpevents(ie).stlas),ie);
	fprintf(afp,' %d\n',ie);
	maxtp = max([twpevents(ie).tp]);
	mintp = min([twpevents(ie).tp]);
	meanamp = mean([twpevents(ie).amp].^0.5);
	for ista=1:length(twpevents(ie).stlas)
		fprintf(cfp,'%s\n',char(twpevents(ie).stnms(ista)));
		fprintf(afp,' 0.0\n');
		dist = twpevents(ie).dists(ista);
		azi = twpevents(ie).azi(ista);
		baz = twpevents(ie).baz(ista);
		stla = twpevents(ie).stlas(ista);
		stlo = twpevents(ie).stlos(ista);
		deg = twpevents(ie).degs(ista);
		amp = twpevents(ie).amp(ista).^0.5;
		stnm = char(twpevents(ie).stnms(ista));
		if length(stnm) == 1
			amp = amp*OBS_amp_bust;
		end
		amp = amp./meanamp;
		tp = twpevents(ie).tp(ista);
%		ph = (maxtp-tp)./periods(ip)/2/pi;
		ph = (tp-mintp)./periods(ip);
		fprintf(afp,' %f  %f  %f  %f  %f  %f\n',dist,azi,baz,deg,stla,stlo);
		fprintf(afp,' %f  %f\n',amp,ph);
	end
end

fprintf(cfp,'%d\n',nfreq);
fprintf(cfp,'%f\n',freq);
fprintf(cfp,'%s\n',detailoutput);
fprintf(cfp,'%s\n',summaroutput);
fprintf(cfp,'%s\n',gridnodes);
fprintf(cfp,'%s\n',ampph_file);
fprintf(cfp,'%s\n',varfile);
fprintf(cfp,'%s\n',mavamp);
fprintf(cfp,'%s\n',fvelarea0);
fprintf(cfp,'%s\n',velarea);
fprintf(cfp,'%d %d %f %f\n',iterlimit,wlambda,dampvel,dampaniso);
fprintf(cfp,'%f\n',unifvel);
fprintf(cfp,'%s\n',kernelfile);

fclose(cfp);
fclose(afp);

end
