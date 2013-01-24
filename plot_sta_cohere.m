% scripts to plot the coherency between OBS and inlands

setup_parameters;

periods = parameters.periods;

badstnms = textread('badsta.lst','%s');

eventcs_path = './CSmeasure/';
csmatfiles = dir(['./CSmeasure/*cs.mat']);

land_land_num = 0;
land_obs_num = 0;
obs_obs_num = 0;
for ie = 1:length(csmatfiles)
	temp = load([eventcs_path,csmatfiles(ie).name]);
	eventcs =  temp.eventcs;
	disp(eventcs.id);
	badstaids = find(ismember(eventcs.stnms,badstnms));
	for ics=1:length(eventcs.CS)
		badstanum = sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids));
		if badstanum == 0
			land_land_num = land_land_num + 1;
			land_land_coheres(:,land_land_num) = eventcs.CS(ics).cohere(:);
		elseif badstanum == 1
			land_obs_num = land_obs_num + 1;
			land_obs_coheres(:,land_obs_num) = eventcs.CS(ics).cohere(:);
		elseif badstanum == 2
			obs_obs_num = obs_obs_num + 1;
			obs_obs_coheres(:,obs_obs_num) = eventcs.CS(ics).cohere(:);
		end
	end
end

figure(17)
clf
N=0;
NP = length(periods);
bins = linspace(0,1,10);
for ip=1:length(periods)
	N = ip;
	subplot(3,length(periods),N)
	hist(land_land_coheres(ip,:),bins)
	title(['land-land, periods: ',num2str(periods(ip))]);
	xlim([0 1])
	N = NP+ip;
	subplot(3,length(periods),N)
	hist(land_obs_coheres(ip,:),bins)
	title(['land-obs, periods: ',num2str(periods(ip))]);
	xlim([0 1])
	N = 2*NP+ip;
	subplot(3,length(periods),N)
	hist(obs_obs_coheres(ip,:),bins)
	title(['obs-obs, periods: ',num2str(periods(ip))]);
	xlim([0 1])
end
