function smap=smoothmap(xi,yi,map,xnode,ynode,D)
% This function is used to make a gaussian smoothing on a grd map

% D is Gaussian average distance

[m n]=size(xi);

midi=floor(m/2);
midj=floor(n/2);

dx=vdist(xi(midi,midj),yi(midi,midj),xi(midi+1,midj),yi(midi+1,midj))/1e3;
dy=vdist(xi(midi,midj),yi(midi,midj),xi(midi,midj+1),yi(midi,midj+1))/1e3;

Nx=floor(D/dx);
Ny=floor(D/dy);

smoothfactor = mean([Nx Ny]);

[smap xi yi] = gridfit(xi(:),yi(:),map(:),xnode,ynode,'smooth',smoothfactor);

end

