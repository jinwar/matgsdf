function smap=smoothmap(xi,yi,map,D)
% This function is used to make a gaussian smoothing on a grd map

% D is average distance

smap = map;
for i=1:length(xi(:))
	dist = distance(xi(i),yi(i),xi,yi);
	dist = deg2km(dist);
	ind = find(dist<D);
	smap(i) = nanmean(map(ind));
end
ind = find(isnan(map(:)));
smap(ind) = NaN;

