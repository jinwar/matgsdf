clear;

setup_parameters
comp = parameters.component;
periods = parameters.periods;

load eikonal_stack_LHZ.mat

for ip=1:length(avgphv)
	filename = ['eikonal_stack_',comp,'_',num2str(periods(ip)),'.xyz'];
	fp = fopen(filename,'w');
	for igrd = 1:length(avgphv(ip).xi(:))
		if ~isnan(avgphv(ip).GV(igrd))
			fprintf(fp,'%f %f %f\n',avgphv(ip).xi(igrd),avgphv(ip).yi(igrd),avgphv(ip).GV(igrd));
		end
	end
end

if ~exist('helmholtz_stack_LHZ.mat','file')
	return;
end

load helmholtz_stack_LHZ.mat;

for ip=1:length(avgphv)
	filename = ['helmholtz_stack_',comp,'_',num2str(periods(ip)),'.xyz'];
	fp = fopen(filename,'w');
	for igrd = 1:length(avgphv(ip).xi(:))
		if ~isnan(avgphv(ip).GV(igrd))
			fprintf(fp,'%f %f %f\n',avgphv(ip).xi(igrd),avgphv(ip).yi(igrd),avgphv(ip).GV_cor(igrd));
		end
	end
end

