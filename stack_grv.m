% Program to stack phase velocity maps from each event

clear;

phase_v_path = './gveikonal/'
r = 0.05;
isfigure = 1;

setup_parameters

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
min_phv_tol = parameters.min_phv_tol;
max_phv_tol = parameters.max_phv_tol;
is_raydense_weight = parameters.is_raydense_weight;
min_event_num = parameters.min_event_num;
err_std_tol = parameters.err_std_tol;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

for ip=1:length(periods)
	avgphv(ip).sumV = zeros(Nx,Ny);
	avgphv(ip).sumweight = zeros(Nx,Ny);
	avgphv(ip).GV_std = zeros(Nx,Ny);
	avgphv(ip).eventnum = zeros(Nx,Ny);
	avgphv(ip).xi = xi;
	avgphv(ip).yi = yi;
	avgphv(ip).xnode = xnode;
	avgphv(ip).ynode = ynode;
	avgphv(ip).period = periods(ip);
end

phvmatfiles = dir([phase_v_path,'/*_',comp,'_gv.mat']);

GV_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
raydense_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));

for ie = 1:length(phvmatfiles)
	temp = load([phase_v_path,phvmatfiles(ie).name]);
	eventphv = temp.eventphv;
	disp(eventphv(1).id);
	evla(ie) = eventphv(ip).evla;
	evlo(ie) = eventphv(ip).evlo;
	for ip=1:length(periods)
        ind = find(eventphv(ip).GV < min_phv_tol);
        eventphv(ip).GV(ind) = NaN;
        ind = find(eventphv(ip).GV > max_phv_tol);
        eventphv(ip).GV(ind) = NaN;
		if eventphv(ip).goodnum./eventphv(ip).badnum < parameters.min_csgoodratio
			disp('too few good cs');
			eventphv(ip).GV(:) = NaN;
        end
        GV_mat(:,:,ie,ip) = eventphv(ip).GV;
		if ~is_raydense_weight
			eventphv(ip).raydense(:) = 1;
		end
        raydense_mat(:,:,ie,ip) = eventphv(ip).raydense;

		ind = find(~isnan(eventphv(ip).GV));
		avgphv(ip).sumV(ind) = avgphv(ip).sumV(ind) + eventphv(ip).GV(ind).*eventphv(ip).raydense(ind);
		avgphv(ip).sumweight(ind) = avgphv(ip).sumweight(ind) + eventphv(ip).raydense(ind);
		avgphv(ip).eventnum(ind) = avgphv(ip).eventnum(ind)+1;
	end
end

for ip=1:length(periods)
	avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).sumweight;
	ind = find(avgphv(ip).eventnum < min_event_num);
	avgphv(ip).GV(ind) = NaN;
end

ori_GV_mat = GV_mat;
% Calculate std, remove the outliers
for ip=1:length(periods)
	for i = 1:Nx
		for j=1:Ny
			avgphv(ip).GV_std(i,j) = nanstd(GV_mat(i,j,:,ip));
			ind = find( abs(GV_mat(i,j,:,ip) - avgphv(ip).GV(i,j)) > err_std_tol*avgphv(ip).GV_std(i,j));
			GV_mat(i,j,ind,ip) = NaN;
		end
	end
end

% calculate the averaged phase velocity again
for ip=1:length(periods)
	avgphv(ip).sumV = zeros(Nx,Ny);
	avgphv(ip).sumweight = zeros(Nx,Ny);
	avgphv(ip).eventnum = zeros(Nx,Ny);
end

for ip = 1:length(periods)
	for ie = 1:length(phvmatfiles)
		GV = GV_mat(:,:,ie,ip);
		raydense = raydense_mat(:,:,ie,ip);
		ind = find(~isnan(GV));
		avgphv(ip).sumV(ind) = avgphv(ip).sumV(ind) + GV(ind).*raydense(ind);
		avgphv(ip).sumweight(ind) = avgphv(ip).sumweight(ind) + raydense(ind);
		avgphv(ip).eventnum(ind) = avgphv(ip).eventnum(ind)+1;
	end
end

for ip=1:length(periods)
	avgphv(ip).GV = avgphv(ip).sumV ./ avgphv(ip).sumweight;
	ind = find(avgphv(ip).eventnum < min_event_num);
	avgphv(ip).GV(ind) = NaN;
end	

save(['eikonal_stack_',comp,'_gv.mat'],'avgphv','GV_mat','ori_GV_mat','evla','evlo');

if isfigure
N=3; M = floor(length(periods)/N)+1;
figure(89)
clf
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	if exist('drawlocal.m','file')
		drawlocal
	end
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
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	if exist('drawlocal.m','file')
		drawlocal
	end
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv(ip).GV_std);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	load seiscmap
	colormap(seiscmap)
	caxis([0 0.4])
end
drawnow;

figure(91)
clf
for ip = 1:length(periods)
	subplot(M,N,ip)
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	if exist('drawlocal.m','file')
		drawlocal
	end
	h1=surfacem(xi,yi,avgphv(ip).sumweight);
	% set(h1,'facecolor','interp');
	title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	meanweight = nanmedian(avgphv(ip).sumweight(:));
	caxis([0 2*meanweight]);
	load seiscmap
	colormap(seiscmap)
end
drawnow;
end
