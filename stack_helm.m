% Program to stack phase velocity maps from each event

clear;

isrecalulatealpha = 0;
fixalpha = 1;
demoip = 4;
isfigure = 1;

phase_v_path = './helmholtz/'
r = 0.05;

setup_parameters

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
min_csgoodratio = parameters.min_csgoodratio;
min_phv_tol = parameters.min_phv_tol;
max_phv_tol = parameters.max_phv_tol;
is_raydense_weight = parameters.is_raydense_weight;
err_std_tol = parameters.err_std_tol;
min_event_num = parameters.min_event_num;
issmoothmap = parameters.issmoothmap;
smooth_wavelength = parameters.smooth_wavelength;


xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

for ip=1:length(periods)
	avgphv(ip).sumV = zeros(Nx,Ny);
	avgphv(ip).sumV_cor = zeros(Nx,Ny);
	avgphv(ip).sumweight = zeros(Nx,Ny);
	avgphv(ip).GV_std = zeros(Nx,Ny);
	avgphv(ip).GV_cor_std = zeros(Nx,Ny);
	avgphv(ip).eventnum = zeros(Nx,Ny);
	avgphv(ip).xi = xi;
	avgphv(ip).yi = yi;
	avgphv(ip).xnode = xnode;
	avgphv(ip).ynode = ynode;
	avgphv(ip).period = periods(ip);
end

phvmatfiles = dir([phase_v_path,'/*_helmholtz_',comp,'.mat']);

GV_cor_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
GV_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
raydense_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));

for ie = 1:length(phvmatfiles)
	temp = load([phase_v_path,phvmatfiles(ie).name]);
	helmholtz = temp.helmholtz;
	if isrecalulatealpha
		for ip=1:length(periods)
			helmholtz(ip).GV_cor = ((helmholtz(ip).GV).^-2 + fixalpha.*helmholtz(ip).amp_term').^-.5;
		end
	end
			
	disp(helmholtz(1).id);
	for ip=1:length(periods)
        ind = find(helmholtz(ip).GV_cor < min_phv_tol);
        helmholtz(ip).GV_cor(ind) = NaN;
        ind = find(helmholtz(ip).GV_cor > max_phv_tol);
        helmholtz(ip).GV_cor(ind) = NaN;
        ind = find(helmholtz(ip).GV < min_phv_tol);
        helmholtz(ip).GV(ind) = NaN;
        ind = find(helmholtz(ip).GV > max_phv_tol);
        helmholtz(ip).GV(ind) = NaN;
		if helmholtz(ip).goodnum./helmholtz(ip).badnum < min_csgoodratio
			disp('not enough good cs measurement');
			helmholtz(ip).GV_cor(:) = NaN;
			helmholtz(ip).GV(:) = NaN;
        end
        GV_cor_mat(:,:,ie,ip) = helmholtz(ip).GV_cor;
        GV_mat(:,:,ie,ip) = helmholtz(ip).GV;
		if ~is_raydense_weight
			helmholtz(ip).raydense(:) = 1;
		end
        raydense_mat(:,:,ie,ip) = helmholtz(ip).raydense;

		ind = find(~isnan(helmholtz(ip).GV_cor));
		avgphv(ip).sumV_cor(ind) = avgphv(ip).sumV_cor(ind) + helmholtz(ip).GV_cor(ind).*helmholtz(ip).raydense(ind);
		avgphv(ip).sumV(ind) = avgphv(ip).sumV(ind) + helmholtz(ip).GV(ind).*helmholtz(ip).raydense(ind);
		avgphv(ip).sumweight(ind) = avgphv(ip).sumweight(ind) + helmholtz(ip).raydense(ind);
		avgphv(ip).eventnum(ind) = avgphv(ip).eventnum(ind)+1;
	end
end

for ip=1:length(periods)
	if ~is_raydense_weight
		avgphv(ip).GV_cor = avgphv(ip).sumV_cor ./ avgphv(ip).eventnum;
		avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).eventnum;
	else
		avgphv(ip).GV_cor = avgphv(ip).sumV_cor ./ avgphv(ip).sumweight;
		avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).sumweight;
	end
	ind = find(avgphv(ip).eventnum < min_event_num);
	avgphv(ip).GV(ind) = NaN;
	avgphv(ip).GV_cor(ind) = NaN;
end


% Calculate std, remove the outliers
for ip=1:length(periods)
	for i = 1:Nx
		for j=1:Ny
			avgphv(ip).GV_std(i,j) = nanstd(GV_mat(i,j,:,ip));
			avgphv(ip).GV_cor_std(i,j) = nanstd(GV_cor_mat(i,j,:,ip));
			ind = find( abs(GV_mat(i,j,:,ip) - avgphv(ip).GV(i,j)) > err_std_tol*avgphv(ip).GV_std(i,j));
			GV_mat(i,j,ind,ip) = NaN;
			ind = find( abs(GV_cor_mat(i,j,:,ip) - avgphv(ip).GV_cor(i,j)) > err_std_tol*avgphv(ip).GV_cor_std(i,j));
			GV_cor_mat(i,j,ind,ip) = NaN;
			avgphv(ip).GV_std(i,j) = nanstd(GV_mat(i,j,:,ip));
			avgphv(ip).GV_cor_std(i,j) = nanstd(GV_cor_mat(i,j,:,ip));
		end
	end
end
% calculate the averaged phase velocity again
for ip=1:length(periods)
	avgphv(ip).sumV = zeros(Nx,Ny);
	avgphv(ip).sumV_cor = zeros(Nx,Ny);
	avgphv(ip).sumweight = zeros(Nx,Ny);
	avgphv(ip).eventnum = zeros(Nx,Ny);
end

for ip = 1:length(periods)
	for ie = 1:length(phvmatfiles)
		raydense = raydense_mat(:,:,ie,ip);
		GV = GV_mat(:,:,ie,ip);
		GV_cor = GV_cor_mat(:,:,ie,ip);
		ind = find(~isnan(GV_cor) & ~isnan(GV));
		if is_raydense_weight
			avgphv(ip).sumV(ind) = avgphv(ip).sumV(ind) + GV(ind).*raydense(ind);
			avgphv(ip).sumV_cor(ind) = avgphv(ip).sumV_cor(ind) + GV_cor(ind).*raydense(ind);
		else
			avgphv(ip).sumV(ind) = avgphv(ip).sumV(ind) + GV(ind);
			avgphv(ip).sumV_cor(ind) = avgphv(ip).sumV_cor(ind) + GV_cor(ind);
		end
		avgphv(ip).sumweight(ind) = avgphv(ip).sumweight(ind) + raydense(ind);
		avgphv(ip).eventnum(ind) = avgphv(ip).eventnum(ind)+1;
	end
end

for ip=1:length(periods)
	if is_raydense_weight
		avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).sumweight;
		avgphv(ip).GV_cor = avgphv(ip).sumV_cor ./ avgphv(ip).sumweight;
	else
		avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).eventnum;
		avgphv(ip).GV_cor = avgphv(ip).sumV_cor ./ avgphv(ip).eventnum;
	end
	ind = find(avgphv(ip).eventnum < min_event_num);
	avgphv(ip).GV(ind) = NaN;
	avgphv(ip).GV_cor(ind) = NaN;
end	

if issmoothmap
	disp(['Smoothing map based on wavelength']);
	for ip=1:length(periods)
		disp(ip);
		D = smooth_wavelength*nanmean(avgphv(ip).GV(:))*periods(ip);
		GV = smoothmap(xi,yi,avgphv(ip).GV,D);
		GV(find(isnan(avgphv(ip).GV))) = NaN;
		avgphv(ip).GV = GV;
		GV_cor = smoothmap(xi,yi,avgphv(ip).GV_cor,D);
		GV_cor(find(isnan(avgphv(ip).GV_cor))) = NaN;
		avgphv(ip).GV_cor = GV_cor;
	end	
end

save(['helmholtz_stack_',comp,'.mat'],'avgphv');

if isfigure
figure(71)
clf
%for ip = 1:length(periods)
ip = demoip
	subplot(2,2,2)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_cor);
	% set(h1,'facecolor','interp');
	title(['stack for corrected phv,','Periods: ',num2str(periods(ip))],'fontsize',15)
	avgv = nanmean(avgphv(ip).GV(:));
	if isnan(avgv)
		continue;
	end
	caxis([avgv*(1-r) avgv*(1+r)])
	colorbar
	load seiscmap
	colormap(seiscmap)
	subplot(2,2,1)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV);
	% set(h1,'facecolor','interp');
	title(['stack for dynamics phv,','Periods: ',num2str(periods(ip))],'fontsize',15)
	avgv = nanmean(avgphv(ip).GV(:));
	if isnan(avgv)
		continue;
	end
	caxis([avgv*(1-r) avgv*(1+r)])
	colorbar
	load seiscmap
	colormap(seiscmap)

	subplot(2,2,3)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_std);
	title(['Original STD,','Periods: ',num2str(periods(ip))],'fontsize',15)
	meanstd = nanmean(avgphv(ip).GV_std(:));
	if ~isnan(meanstd)
		caxis([0 2*meanstd])
	end
	colorbar
	load seiscmap
	colormap(seiscmap)

	subplot(2,2,4)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_cor_std);
	title(['corrected STD,','Periods: ',num2str(periods(ip))],'fontsize',15)
	if ~isnan(meanstd)
		caxis([0 2*meanstd])
	end
	colorbar
	load seiscmap
	colormap(seiscmap)

N=3; M = floor(length(periods)/N)+1;
figure(89)
clf
title('stack for dynamics phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	avgv = nanmean(avgphv(ip).GV(:));
	if isnan(avgv)
		continue;
	end
	caxis([avgv*(1-r) avgv*(1+r)])
	colorbar
	load seiscmap
	colormap(seiscmap)
end
drawnow;

figure(90)
clf
title('Std for dynamics phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_std);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
	meanstd = nanmean(avgphv(ip).GV_std(:));
	if ~isnan(meanstd)
		caxis([0 2*meanstd])
	end
end
drawnow;

figure(91)
clf
title('stack for structure phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_cor);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	avgv = nanmean(avgphv(ip).GV(:));
	if isnan(avgv)
		continue;
	end
	caxis([avgv*(1-r) avgv*(1+r)])
	colorbar
	load seiscmap
	colormap(seiscmap)
end
drawnow;

figure(92)
clf
title('Std for structure phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_cor_std);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
	meanstd = nanmean(avgphv(ip).GV_std(:));
	if ~isnan(meanstd)
		caxis([0 2*meanstd])
	end
%	caxis([0 0.5])
end
drawnow;

figure(93)
clf
title('diff phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_cor-avgphv(ip).GV);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
%	caxis([0 0.5])
end
drawnow;

figure(93)
clf
title('diff phv')
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_cor_std-avgphv(ip).GV_std);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
%	caxis([0 0.5])
end
drawnow;

figure(95)
clf
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).sumweight);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
end
drawnow;
end
