clear;

period = 40;
control_file = ['all.0',num2str(period)];
phamp_file = ['phampcor.0',num2str(period),'.inp'];

cfp = fopen(control_file,'r');
stemp = fgetl(cfp);
eventnum = sscanf(stemp,'%d');

for ie = 1:eventnum
	stemp = fgetl(cfp);
	dtemp = sscanf(stemp,'%d %d');
	stanum(ie) = dtemp(1);
	for ista = 1:stanum(ie)
		stemp = fgetl(cfp);
	end
end	
fclose(cfp);

dfp = fopen(phamp_file,'r');
for ie = 1:eventnum
	stemp = fgetl(dfp);
	eventid = sscanf(stemp,'%d');
	if eventid ~=ie
		disp('something wrong!');
	end
	event(ie).id = eventid;
	for ista = 1:stanum(ie)
		stemp = fgetl(dfp);
		event(ie).sta(ista).bgtime = sscanf(stemp,'%f');
		stemp = fgetl(dfp);
		dtemp = sscanf(stemp,'%f %f %f %f %f %f');
		event(ie).sta(ista).dist = dtemp(1);
		event(ie).sta(ista).azi = dtemp(2);
		event(ie).sta(ista).baz = dtemp(3);
		event(ie).sta(ista).deg = dtemp(4);
		event(ie).sta(ista).lat = dtemp(5);
		event(ie).sta(ista).lon = dtemp(6);
		stemp = fgetl(dfp);
		dtemp = sscanf(stemp,'%f %f');
		event(ie).sta(ista).amp = dtemp(1);
		event(ie).sta(ista).ph = dtemp(2);
	end
end

% normalize the phase

alldist = [];
allphs = [];
for ie = 1:eventnum
	dists = [event(ie).sta.dist];
	phs = [event(ie).sta.ph];
	[mindist mindistid] = min(dists);
	dphs = phs - phs(mindistid);
	ddists = dists - dists(mindistid);
	alldist = [alldist;ddists(:)];
	allphs = [allphs;dphs(:)];
end
		
figure(23)
clf
hold on
plot(alldist,allphs,'x');
[x y] = ginput(2);
plot(x,y,'r');

phv = abs(diff(x))/abs(diff(y)*period)
