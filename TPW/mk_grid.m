% Script to print out grid nodes needed for fortran two plane wave method
% NJA, LDEO 2013
%
% Format is as follows and moves from left to right
% ncol x nrow
% # of grid nodes
% latitude longitude
% corner points x4
% ncol
% col_spacing row_spacing (in km)

% Parameters
GRID_file = 'PNG_gridnodes';
lalim = [-11.2 -7.8];
lolim = [148.8 151.5];
gridsize = 0.1;
gridkm = gridsize*111.32;
gridla = lalim(2):-gridsize:lalim(1);
gridlo = lolim(1):gridsize:lolim(2);

ncol = length(gridlo);
nrow = length(gridla);
ngrid = ncol*nrow;

head = sprintf('grid %sx%s',num2str(ncol),num2str(nrow));

fid = fopen(GRID_file,'w');
if (fid == -1)
    error (['    Cannot open file: ', F_OUT]);
end



% Print header information
fprintf(fid,'%s\n',head);
fprintf(fid,'%d\n',ngrid);

% Print out grid nodes
for ilo = (1:ncol)
    for ila = (1:nrow)
        fprintf(fid,'%4.1f %4.1f\n',[gridla(ila) gridlo(ilo)]);
    end
end

% Print out corner points
fprintf(fid,'%4.1f %4.1f\n',[lalim(1) lolim(1)]);
fprintf(fid,'%4.1f %4.1f\n',[lalim(2) lolim(1)]);
fprintf(fid,'%4.1f %4.1f\n',[lalim(2) lolim(2)]);
fprintf(fid,'%4.1f %4.1f\n',[lalim(1) lolim(2)]);


% print out ncol
fprintf(fid,'%2i\n',ncol);

% print out col_spacing krow_spacing (in km)
fprintf(fid,'%4.1f %4.1f\n',[gridkm gridkm]);
fclose(fid);
disp('Done Done')
