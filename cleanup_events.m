% script to remove events that are close in the time that may interference with each other.
% written by Ge Jin,jinwar@gmail.com


clear;

eventmatpath = './eventmat/';
csmatpath = './CSmeasure/';
eikonalpath = './eikonal/';

% Setup parameters
setup_parameters

% Setup Error Codes for Bad data
setup_ErrorCode

lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
gridsize = gridsize*3;
xnode = [lalim(1),mean(lalim),lalim(2)];
ynode = [lolim(1),mean(lolim),lolim(2)];
[xi yi] = meshgrid(xnode,ynode);

matfiles = dir([eventmatpath,'/*_',parameters.component,'.mat']);
for ie = 1:length(matfiles)
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;
	disp(event.id)
	evotimes(ie) = event.otime;
	evlas(ie) = event.evla;
	evlos(ie) = event.evlo;
	evids(ie) = {event.id};
	isgood(ie) = 1;
	dist = deg2km(distance(xi,yi,evlas(ie),evlos(ie)));
	win_start(ie) = evotimes(ie) + min(dist(:)/5);
	win_end(ie) = evotimes(ie) + max(dist(:)/2);
	if ~isfield(event,'stadata')
		isgood(ie)=0;
	elseif length(event.stadata)<5
		isgood(ie)=0;
	end
end % end of event loop

for ie = 1:length(evlas)
	for je = 1:length(evlas)
		if ie == je
			continue;
		end
		if win_start(ie) > win_start(je) && win_start(ie) < win_end(je)
			isgood(ie) = 0;
			isgood(je) = 0;
		end
		if win_end(ie) > win_start(je) && win_end(ie) < win_end(je)
			isgood(ie) = 0;
			isgood(je) = 0;
		end
	end
end

disp('Bad events:')
badind = find(isgood == 0);
for ie = badind
	disp(evids(ie));
end

%com = input('Do you want to delete these events? y/n','s');
com = 'y';

if com == 'y'
for ie = badind
	delete([eventmatpath,char(evids(ie)),'*.mat'])
	delete([csmatpath ,char(evids(ie)),'*.mat'])
	delete([eikonalpath ,char(evids(ie)),'*.mat'])
end
end

