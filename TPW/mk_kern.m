% Script to run seizmo script makekernels to make amplitude and phase
% sensitivity kernels for Yang and Forsyth codes.
% NJA, Summer 2013
clear all

% Parameters
periods = [20 25 32 40 50 60];
refvs = [3.53 3.61 3.7 3.76 3.82 3.85];

for ip = 1:6
	period = periods(ip);
	smlength = 40;


	range = [1/period 1/period];
	width = 0.10;
	offset = 0.1;

	% First, get some frequency bands:
	band=filter_bank(range,'variable',width,offset);

	% Second, get some phase velocities at those frequencies:
	%phvel=prem_dispersion(band(:,1));
	phvel=refvs(ip);

	% Third, use the beat length as a window width:
	L_beat=1./(band(:,3)-band(:,2));

	% Fourth, modify that window width by an empirically derived
	% power law to get a more typical window width:
	win=L_beat'*(2.5+1000*band(:,1).^2);

	% Now get the kernels for a windowed sinusoid sampled at 1Hz, tapered
	% on the outer 20% of the window, and zero-padded to a length of
	% 9000s starting at -2000s.  The kernels are sampled in a 5000km grid
	% centered on the receiver, grid points are every 10km and a Gaussian
	% smoother with a characteristic distance of 100km is applied.  The
	% kernel names are appended with '.example'.
	 makekernels(band(:,1),phvel,1,win*[0 1],0.2,[-2000 7000],1000,10,...
		smlength,[],'.example');
 end
