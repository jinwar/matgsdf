% Read in the eventcs structures and apply eikonal tomography on each event.
% Written by Ge Jin, jinwar@gmail.com
% 2013.1.16
%
clear

% debug setting
isfigure = 0;
isdisp = 0;
is_overwrite = 0;

% input path
eventcs_path = './CSmeasure/';
% output path
eikonl_output_path = './gveikonal/';
if ~exist(eikonl_output_path)
	mkdir(eikonl_output_path);
end

% setup parameters
setup_parameters

comp = parameters.component;
lalim=parameters.lalim;
lolim=parameters.lolim;
gridsize=parameters.gridsize;
periods = parameters.periods;
raydensetol=parameters.raydensetol;
smweight0 = parameters.smweight0;
Tdumpweight0 = parameters.Tdumpweight;
Rdumpweight0 = parameters.Rdumpweight;
fiterrtol = parameters.fiterrtol;
dterrtol = parameters.dterrtol;
isRsmooth = parameters.isRsmooth;
inverse_err_tol = parameters.inverse_err_tol;
min_amp_tol  = parameters.min_amp_tol;

% setup useful variables
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

% Setup universal smoothing kernel
disp('initial the smoothing kernel')
tic
    [i,j] = ndgrid(1:Nx,2:(Ny-1));
    ind = j(:) + Ny*(i(:)-1);
    dy = diff(ynode);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
                    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny,Nx*Ny);

    [i,j] = ndgrid(2:(Nx-1),1:Ny);
    ind = j(:) + Ny*(i(:)-1);
    dx = diff(xnode);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
            [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny,Nx*Ny)];

    F=sparse(Nx*Ny*2*2,Nx*Ny*2);
    for n=1:size(Areg,1)
        ind=find(Areg(n,:)~=0);
        F(2*n-1,2*ind-1)=Areg(n,ind);
        F(2*n,2*ind)=Areg(n,ind);
    end
toc

% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	disp('Found Bad stations:')
	disp(badstnms)
end

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
for ie = 1:length(csmatfiles)
%for ie = 30
	clear eventphv 
	% read in data and set up useful variables
	temp = load([eventcs_path,csmatfiles(ie).name]);
	eventcs =  temp.eventcs;
	disp(eventcs.id)
	evla = eventcs.evla;
	evlo = eventcs.evlo;

	matfilename = [eikonl_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat'];
	if exist(matfilename,'file') && ~is_overwrite
		disp(['Exist ',matfilename,', skip']);
		continue;
	end


	if exist('badstnms','var')
		badstaids = find(ismember(eventcs.stnms,badstnms));
	else
		badstaids = [];
	end

	% Build the rotation matrix
	razi = azimuth(xi+gridsize/2,yi+gridsize/2,evla,evlo)+180;
	R = sparse(2*Nx*Ny,2*Nx*Ny);
	for i=1:Nx
		for j=1:Ny
			n=Ny*(i-1)+j;
			theta = razi(i,j);
			R(2*n-1,2*n-1) = cosd(theta);
			R(2*n-1,2*n) = sind(theta);
			R(2*n,2*n-1) = -sind(theta);
			R(2*n,2*n) = cosd(theta);
		end
	end

	clear rays
	for ics = 1:length(eventcs.CS)
		rays(ics,1) = eventcs.stlas(eventcs.CS(ics).sta1);
		rays(ics,2) = eventcs.stlos(eventcs.CS(ics).sta1);
		rays(ics,3) = eventcs.stlas(eventcs.CS(ics).sta2);
		rays(ics,4) = eventcs.stlos(eventcs.CS(ics).sta2);
	end

	% Build the kernel
	disp('Buildling up ray path kernel')
	tic
		mat=kernel_build(rays,xnode,ynode);
	toc

	% build dumping matrix for ST
	dumpmatT = R(2:2:2*Nx*Ny,:);
	
	% build dumping matrix for SR
	dumpmatR = R(1:2:2*Nx*Ny-1,:);
	
	% Loop through the periods
	for ip = 1:length(periods)
		dt = zeros(length(eventcs.CS),1);
		w = zeros(length(eventcs.CS),1);
		for ics = 1:length(eventcs.CS)
			if eventcs.CS(ics).isgood(ip) > 0 
				dt(ics) = eventcs.CS(ics).dtg(ip);
				w(ics) = 1;
			else
				dt(ics) = 0;
				w(ics) = 0;
			end
			if sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids)) > 0
				w(ics) = 0;
			end
		end
		W = sparse(diag(w));

		% Normalize smoothing kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;

		% Normalize dumping matrix for ST
		NR=norm(dumpmatT,1);
		NA=norm(W*mat,1);
		dumpweightT = Tdumpweight0*NA/NR;
		
		% Normalize dumping matrix for SR
		NR=norm(dumpmatR,1);
		NA=norm(W*mat,1);
		dumpweightR = Rdumpweight0*NA/NR;

		% Set up matrix on both side
		if isRsmooth
            A=[W*mat;smweight*F*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
        else
            A=[W*mat;smweight*F;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
        end

		avgv = eventcs.avgphv(ip);
        rhs=[W*dt;zeros(size(F,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
        
		% Least square inversion
        phaseg=(A'*A)\(A'*rhs);
	        
        % Iteratively down weight the measurement with high error
		niter=0;
		ind = find(diag(W)==0);
		if isdisp
			disp(['Before iteration'])
			disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
			disp(['Bad Measurement Number: ', num2str(length(ind))]);
		end
        niter=1;
        while niter < 2
            niter=niter+1;
            err = mat*phaseg - dt;
			err = W*err;
            %            err = W*err;
            stderr=std(err);
            if stderr > dterrtol
                stderr = dterrtol;
            end
            for i=1:length(err)
                if abs(err(i)) > inverse_err_tol*stderr
                    W(i,i)=0;
                end
            end
            ind = find(diag(W)==0);
			if isdisp
				disp('After:')
				disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
				disp(['Bad Measurement Number: ', num2str(length(ind))]);
			end
            
            % Rescale the smooth kernel
            NR=norm(F,1);
            NA=norm(W*mat,1);
            smweight = smweight0*NA/NR;
            
            % rescale dumping matrix for St
            NR=norm(dumpmatT,1);
            NA=norm(W*mat,1);
            dumpweightT = Tdumpweight0*NA/NR;
            
            % rescale dumping matrix for SR
            NR=norm(dumpmatR,1);
            NA=norm(W*mat,1);
            dumpweightR = Rdumpweight0*NA/NR;
            
            if isRsmooth
                A=[W*mat;smweight*F*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
            else
                A=[W*mat;smweight*F;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
            end
            rhs=[W*dt;zeros(size(F,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
            phaseg=(A'*A)\(A'*rhs);
        end	

        % Calculate the kernel density
        %sumG=sum(abs(mat),1);
        ind=1:Nx*Ny;
        rayW = W;
        rayW(find(rayW>1))=1;
        raymat = rayW*mat;
        sumG(ind)=sum((raymat(:,2*ind).^2+raymat(:,2*ind-1).^2).^.5,1);
        clear raydense
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                raydense(i,j)=sumG(n);
            end
        end
        
        %        disp(' Get rid of uncertainty area');
        fullphaseg = phaseg;
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                if raydense(i,j) < raydensetol %&& ~issyntest
                    phaseg(2*n-1)=NaN;
                    phaseg(2*n)=NaN;
                end
            end
        end

		% Change phaseg into phase velocity
		for i=1:Nx
			for j=1:Ny
				n=Ny*(i-1)+j;
				GVx(i,j)= phaseg(2*n-1);
				GVy(i,j)= phaseg(2*n);
			end
		end
		GV=(GVx.^2+GVy.^2).^-.5;
		% Get rid of uncertain area

		% save the result in the structure
		eventphv(ip).rays = rays;
		eventphv(ip).w = diag(W);
		eventphv(ip).goodnum = length(find(w>0));
		eventphv(ip).badnum = length(find(w==0));
		eventphv(ip).dt = dt;
		eventphv(ip).GV = GV;
		eventphv(ip).GVx = GVx;
		eventphv(ip).GVy = GVy;
		eventphv(ip).raydense = raydense;
		eventphv(ip).lalim = lalim;
		eventphv(ip).lolim = lolim;
		eventphv(ip).gridsize = gridsize;
		eventphv(ip).id = eventcs.id;
		eventphv(ip).evla = eventcs.evla;
		eventphv(ip).evlo = eventcs.evlo;
		eventphv(ip).period = periods(ip);
		disp(['Period:',num2str(periods(ip)),', Goodnum:',num2str(eventphv(ip).goodnum),...
				'Badnum:',num2str(eventphv(ip).badnum)]);
	end % end of periods loop
	if isfigure
		N=3; M = floor(length(periods)/N) +1;
		figure(88)
		clf
		for ip = 1:length(periods)
			subplot(M,N,ip)
			ax = worldmap(lalim, lolim);
			set(ax, 'Visible', 'off')
			h1=surfacem(xi,yi,eventphv(ip).GV);
			% set(h1,'facecolor','interp');
%			load pngcoastline
%			geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)
			title(['Periods: ',num2str(periods(ip))],'fontsize',15)
			avgv = nanmean(eventphv(ip).GV(:));
			if isnan(avgv)
				continue;
			end
			r = 0.1;
			caxis([avgv*(1-r) avgv*(1+r)])
			colorbar
			load seiscmap
			colormap(seiscmap)
		end
		drawnow;
	end
	matfilename = [eikonl_output_path,'/',eventcs.id,'_eikonal_',comp,'_gv.mat'];
	save(matfilename,'eventphv');
	disp(['Save the result to: ',matfilename])
end % end of loop ie
