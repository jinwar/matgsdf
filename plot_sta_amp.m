% code to plot the average amplitude of all stations

clear
setup_parameters;
comp = parameters.component;

CSfiles = dir(['CSmeasure/*_',comp,'*.mat']);

% Gather information
stnms = {};
stainfo = [];
for ie = 1:length(CSfiles)
	clear eventcs amps
	load(fullfile('CSmeasure',CSfiles(ie).name));
	disp(CSfiles(ie).name)
	for ista = 1:length(eventcs.stnms)
		amps(ista,:) = sqrt(eventcs.autocor(ista).amp);
	end
	meanamp = mean(amps,1);

	for ista = 1:length(eventcs.stnms)
		if ismember(eventcs.stnms(ista),stnms)
			staid = find(ismember(stnms,eventcs.stnms(ista)));
			stainfo(staid).norm_amp(end+1,:) = amps(ista,:)./meanamp;
			stainfo(staid).ori_amp(end+1,:) = amps(ista,:);
		else
			stnms(end+1) = eventcs.stnms(ista);
			stainfo(end+1).stnm = eventcs.stnms(ista);
			stainfo(end).norm_amp(1,:) = amps(ista,:)./meanamp;
			stainfo(end).ori_amp(1,:) = amps(ista,:);
		end
	end
end

save('png_amp_stainfo.mat','stainfo');

load png_amp_stainfo

avg_band = 3:6;

badstnms = textread('badampsta.lst','%s');
OBSstnms = textread('OBS.lst','%s');

% calculate the means
for ista = 1:length(stainfo)
	stainfo(ista).meanamp = mean(stainfo(ista).norm_amp,1);
	stainfo(ista).avgmean = mean(stainfo(ista).meanamp(avg_band));
	if ismember(stainfo(ista).stnm,badstnms)
		stainfo(ista).isgood = 0;
	else
		stainfo(ista).isgood = 1;
	end
	if ismember(stainfo(ista).stnm,OBSstnms)
		stainfo(ista).isOBS = 1;
	else
		stainfo(ista).isOBS = 0;
	end
end

%calculate the correction
isgood = [stainfo.isgood];
isOBS = [stainfo.isOBS];
onland_ind = find(isgood==1 & isOBS == 0);
onland_avg = mean([stainfo(onland_ind).avgmean]);
OBS_ind = find(isgood==1 & isOBS ==1);
OBS_avg = mean([stainfo(OBS_ind).avgmean]);
OBS_cor = onland_avg/OBS_avg;

%% making plots
% Before the correction
figure(23)
clf
subplot(2,1,1)
hold on
for ip = avg_band
	for ista = 1:length(stainfo)
		norm_amp = stainfo(ista).norm_amp;
		x = ones(size(norm_amp,1),1)*ista;
		plot(x,norm_amp(:,ip),'x');
		errorbar(ista, mean(norm_amp(:,ip)), std(norm_amp(:,ip)),'ro','markerfacecolor','r');
		if stainfo(ista).isgood
			plot(ista, stainfo(ista).avgmean,'ro','markerfacecolor','g');
		end
		maxy = 3;
		ylim([0 maxy]);
%		text(ista,-maxy/10,char(stainfo(ista).stnm),'rotation',90);
		text(ista,-maxy/50,char(stainfo(ista).stnm),'rotation',-90);
		set(gca, 'XTickLabel','')
	end
end
subplot(2,1,2)
isgood = [stainfo.isgood];
ind = find(isgood);
hist([stainfo(ind).avgmean],10);
xlim([0 2])

% Gather information
stnms = {};
stainfo = [];
for ie = 1:length(CSfiles)
	clear eventcs amps
	load(fullfile('CSmeasure',CSfiles(ie).name));
	disp(CSfiles(ie).name)
	for ista = 1:length(eventcs.stnms)
		if ismember(eventcs.stnms(ista),OBSstnms)
			amps(ista,:) = eventcs.autocor(ista).amp.^.5*OBS_cor;
		else
			amps(ista,:) = eventcs.autocor(ista).amp.^.5;
		end
	end
	meanamp = mean(amps,1);

	for ista = 1:length(eventcs.stnms)
		if ismember(eventcs.stnms(ista),stnms)
			staid = find(ismember(stnms,eventcs.stnms(ista)));
			stainfo(staid).norm_amp(end+1,:) = amps(ista,:)./meanamp;
			stainfo(staid).ori_amp(end+1,:) = amps(ista,:);
		else
			stnms(end+1) = eventcs.stnms(ista);
			stainfo(end+1).stnm = eventcs.stnms(ista);
			stainfo(end).norm_amp(1,:) = amps(ista,:)./meanamp;
			stainfo(end).ori_amp(1,:) = amps(ista,:);
		end
	end
end

% calculate the means
for ista = 1:length(stainfo)
	stainfo(ista).meanamp = mean(stainfo(ista).norm_amp,1);
	stainfo(ista).avgmean = mean(stainfo(ista).meanamp(avg_band));
	if ismember(stainfo(ista).stnm,badstnms)
		stainfo(ista).isgood = 0;
	else
		stainfo(ista).isgood = 1;
	end
	if ismember(stainfo(ista).stnm,OBSstnms)
		stainfo(ista).isOBS = 1;
	else
		stainfo(ista).isOBS = 0;
	end
end

save('png_amp_stainfo_cor.mat','stainfo');

clear stainfo

load('png_amp_stainfo_cor');

% after the correction
figure(24)
clf
hold on
title('After Correction');
subplot(2,1,1)
hold on
for ip = avg_band
	for ista = 1:length(stainfo)
		norm_amp = stainfo(ista).norm_amp;
		x = ones(size(norm_amp,1),1)*ista;
		plot(x,norm_amp(:,ip),'x');
		errorbar(ista, mean(norm_amp(:,ip)), std(norm_amp(:,ip)),'ro','markerfacecolor','r');
		if stainfo(ista).isgood
			plot(ista, stainfo(ista).avgmean,'ro','markerfacecolor','g');
		end
		maxy = 3;
		ylim([0 maxy]);
%		text(ista,-maxy/10,char(stainfo(ista).stnm),'rotation',90);
		text(ista,-maxy/50,char(stainfo(ista).stnm),'rotation',-90);
		set(gca, 'XTickLabel','')
	end
end
subplot(2,1,2)
isgood = [stainfo.isgood];
ind = find(isgood);
hist([stainfo(ind).avgmean],10);
xlim([0 2])
