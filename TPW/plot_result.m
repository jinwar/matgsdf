clear;

setup_parameters;
periods = parameters.periods;

for ip = 1:6
period = periods(ip);
nodefile = 'PNG_gridnodes';
phasefile = ['temp.0',num2str(period),'.inp.sa360kern'];
gridsize = 0.1;
r = 0.05;

nodefp = fopen(nodefile,'r');
temp = fgetl(nodefp);
nodenum = sscanf(fgetl(nodefp),'%d');
for i=1:nodenum
	ntemp = sscanf(fgetl(nodefp),'%f %f\n');
	lat(i) = ntemp(1);
	lon(i) = ntemp(2);
end
fclose(nodefp);
lalim = [min(lat) max(lat)];
lolim = [min(lon) max(lon)];
latgrid = lalim(1):gridsize:lalim(2);
longrid = lolim(1):gridsize:lolim(2);

%phvfp = fopen(phasefile,'r');
%
%for i=1:283
%	stemp = fgetl(phvfp);
%end
data = load(phasefile);

%for i=1:nodenum
%	stemp = fgetl(phvfp);
%	ntemp = sscanf(stemp,'%d %f %f');
%	phv(i) = ntemp(2);
%	phvvar(i) = ntemp(3);
%end

phv = data(:,2);
phvvar = data(:,3);

[phv_surf xi yi] = gridfit(lat,lon,phv,latgrid,longrid);
[phvvar_surf xi yi] = gridfit(lat,lon,phvvar,latgrid,longrid);

%figure(32+ip)
%clf
%ax = worldmap(lalim,lolim);
%set(ax, 'Visible', 'off')
%h1=surfacem(xi,yi,phv_surf);
%drawpng
%caxis(mean(phv_surf(:))*[1-r 1+r]);
%colorbar
%load seiscmap
%colormap(seiscmap);
%title(num2str(period),'fontsize',20);

figure(31)
subplot(2,3,ip)
ax = worldmap(lalim,lolim);
set(ax, 'Visible', 'off')
h1=surfacem(xi,yi,phv_surf);
drawpng
caxis(mean(phv_surf(:))*[1-r 1+r]);
colorbar
load seiscmap
colormap(seiscmap);
title(num2str(period),'fontsize',20);

end

%figure(33)
%clf
%ax = worldmap(lalim,lolim);
%set(ax, 'Visible', 'off')
%h1=surfacem(xi,yi,phvvar_surf);
%drawpng
%colorbar
%load seiscmap
%colormap(seiscmap);
%caxis([0 median(phvvar_surf(:))*2]);
