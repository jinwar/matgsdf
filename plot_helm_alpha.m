% Program to plot the best alpha distribution for all frequency bands.

clear;

setup_parameters;
phase_v_path = './helmholtz/'
comp = parameters.component;
periods = parameters.periods;

phvmatfiles = dir([phase_v_path,'/*_helmholtz_',comp,'.mat']);

for ie = 1:length(phvmatfiles)
	disp(phvmatfiles(ie).name);
	temp = load([phase_v_path,phvmatfiles(ie).name]);
	helmholtz = temp.helmholtz;
	alpha(:,ie) = [helmholtz.bestalpha]';
end

figure(12)
clf
N=3; M = floor(length(periods)/N)+1;
for ip=1:length(periods);
	subplot(M,N,ip)
	hist(alpha(ip,:));
	title(num2str(periods(ip)));
end
	
