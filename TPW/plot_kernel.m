clear;

kernel_filename = 'kernel.032s040km3954v.example';

temp = textread(kernel_filename);
data = temp(3:end,:);

x = data(:,1);
y = data(:,2);
ph = data(:,3);
amp = data(:,4);

xnode = min(x):10:max(x);
ynode = min(y):10:max(y);

[ph_surf,xi,yi] = gridfit(x,y,ph,xnode,ynode);
[amp_surf,xi,yi] = gridfit(x,y,amp,xnode,ynode);

figure(44)
clf
surface(xi,yi,ph_surf);
shading flat;

figure(45)
clf
surface(xi,yi,amp_surf);
shading flat;
