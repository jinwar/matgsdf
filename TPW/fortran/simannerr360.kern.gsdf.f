c  simannerr360.kern.f
c  invert for 2-D phase velocities and  station corrections

c  include the sensitivtity kernel for both phase and amplitude
c  the sensitivity kernel for phase is used to calcualte the phase at each  station
c  the sensitivity kernel is calculated forehand. With the sensitivity 
c  kernel data, we can know the sensitivity at each node by interpolating the 
c  sensitivity kernel. 


c  this version includes annealing step.  Inside annealing process, two plane wave interference 
c (the previous mehtod) is used to 
c  calculated the predicted phase and amplitude for each station with 
c six wave parameters and velocity structure, since the annealing with kernel is too slow 

c   compile the code using "gfortran simannerr360.kern.f anneal.f ran1.f -o simannerr360.kern"

       parameter (maxnfreq=1, maxnsta=120, maxnstatype = 1500,
     1    maxpts=20000, nparam = 5000,
     1    maxnobs = 50000, maxnodes = 2000, maxevnts = 500,
     1    maxnxints = 501, maxnyints = 501,ndeg = 81 )

cxx  variables for data 

      real*4 dph(maxevnts,maxnsta)
      real*4 staph(maxevnts,maxnsta,maxnfreq)
      real*4 staamp(maxevnts,maxnsta,maxnfreq)
      real*4 amprms(maxevnts,maxnfreq)
      real*4 streal(maxevnts,maxnsta,maxnfreq)
      real*4 stimag(maxevnts,maxnsta,maxnfreq)
      real*4 bazi(maxevnts,maxnsta), stadelt(maxevnts,maxnsta)
      real*4 d(maxnobs)
      real*4 beg(maxevnts,maxnsta)
      real*4 stadist(maxevnts,maxnsta), staazi(maxevnts,maxnsta)
      real*4 tazim, delta, tazimref, bazimnd(maxnodes)

cxx  variables for the inversion
      real*4 freq(maxnfreq)

c      real*4 g(maxnobs,nparam)
      real*4,allocatable::g(:,:)
      real*4 stddevdata(maxevnts)

      double precision change(nparam), gtdcmm(nparam),ddd
c      double precision gtg(nparam,nparam),gtginv(nparam,nparam)

      double precision, allocatable::gtg(:,:)
      double precision, allocatable::gtginv(:,:)

cc      double precision savegtg(nparam,nparam)

cc      real*4 gtg(nparam,nparam), gtginv(nparam,nparam),ddd
cc      real*4 savegtg(nparam,nparam)

      real*8 chmax
      real*4  gtd(nparam), stddev(nparam), covinv(nparam)
      real*4 origmod(nparam),crrntmod(nparam)


cxx  variables for the coordinates 

      real*4 ysta(maxevnts,maxnsta), xsta(maxevnts,maxnsta)
      real*4 rloc(maxevnts,maxnsta), azloc(maxevnts,maxnsta)
      real*4 attnfac(maxnfreq), xmin(maxevnts,ndeg)
      real*4 boxlat(4), boxlon(4),applat,applon
      real*4 nodelat(maxnodes),nodelon(maxnodes), nodevel(maxnodes)
      real*4 nodecos2(maxnodes), nodesin2(maxnodes)
      real*4 xbox(maxevnts,4),ybox(maxevnts,4)
      real*4 xnode(maxevnts,maxnodes),ynode(maxevnts,maxnodes)
      real*4 stalat(maxevnts,maxnsta),stalon(maxevnts,maxnsta)

cxx  variables for the kernels
      real*4 adistsq(maxnodes),wgttemp(maxnodes)

       real*4,allocatable::wgtnode1(:,:,:),ampwgtnode1(:,:,:)

c      real*4 ampwgtnode1(maxnsta,maxnodes,ndeg)

      real*4 sensitivity(maxnxints,maxnyints)
      real*4 ampsens(maxnxints,maxnyints)
      real*4 dxnode,dynode
      real*4   xbegkern,dxkern
      real*4   ybegkern,dykern
      integer*4 nykern
      integer*4 nxkern

cxx variables in the inversion 

      real*4 dtime(maxnsta),avslow(maxnsta)     
      real*4 dtime1(maxnsta),dtime2(maxnsta),avslow1(maxnsta),
     1         avslow2(maxnsta)
      real*4 appvel(maxnodes),vage(maxnodes)
      real*4 cos2node(maxevnts,maxnodes),sin2node(maxevnts,maxnodes)
      real*4 startamp1(maxevnts),startamp2(maxevnts),stazim1(maxevnts)
      real*4 stazim2(maxevnts),stphase1(maxevnts),stphase2(maxevnts)
       real*4 pv(6),p(7,6),smsft(7),pb(6),annlscl(6),ppb(6),pppb(6)
      real*4 minstd, rmsphase(maxevnts),rmsamp(maxevnts)
      real*4 sortrms(maxevnts)
      real*4 rmsdata(maxevnts)
      real*4  damp1per(maxevnts,maxnsta), damp2per(maxevnts,maxnsta)
      real*4  dphase1(maxevnts,maxnsta),dphase2(maxevnts,maxnsta)
      real*4   phase1(maxevnts,maxnsta), phase2(maxevnts,maxnsta)
      real*4 phase(maxnsta,ndeg),dphase(maxnsta,ndeg),
     1      dampper(maxnsta,ndeg)
      real*4 unifvel

cxx  variables for station corrections
      real*4 ampmult(maxnstatype)
      real*4   phcor(maxnstatype)
      integer*4 istatype(maxnsta),nevntsta(maxnsta)
      integer*4 itypeflag(maxnstatype),iptype(maxnstatype)
      integer*4 istaflag(maxnsta)
      integer*4  sttypename(maxnstatype)
      integer*4  istanum(maxevnts,maxnsta)
      integer*4 nsta(maxevnts), iref(maxevnts)
      integer*4 nfreq,nstapts(maxnsta),indx(nparam),nstacor,nobs
      integer*4 istacor(maxnsta), nnodes, nevents
      integer*4  idnum(maxevnts)
      integer*4  istavar(maxnsta)
      character*6 idstaname

cx  varibles for file names and misc
      integer blank, iblank1, iblank2, iblank3, iblank4
      character*70 foutput0, fsummary0, fvariance0
      character*100 foutput, fn(maxevnts,maxnsta) , fsummary
      character dirpath *34, evtdatatime *14
      character fnname(maxevnts,maxnsta) *70
      character*70 finvrsnodes, fftinput,fvariance,fresmax
      character*70 fmaxavamp,ftemp,fvelarea,fvelarea0
      character*100 sensfn,gridvelfn 
      character*100 dummy 


      common /residua/ d,rloc,azloc,freq,xsta,dtime,avslow,
     1     streal,stimag,stddevdata,iref,nsta,iev,naddat,ifreq,
     1  xbox,ybox,ysta,nxkern,nykern,dxkern,dykern,dxnode,dynode,
     1 sensitivity,ampsens,unifvel,appvel,ampmult,istavar,istanum,
     1  xnode, ynode,nnodes,xmin,
     1   phase,dphase,dampper,phcor,ntype,iptype,istatype

c     common /gmatrix/  gtg,gtginv

      real misfit
      external misfit    
      integer findid
      integer findstatype
  
      pi = 3.1415928
      convdeg = 3.1415928/180.
      circ = 6371.*3.1415928/180.
      twopi = 3.1415928*2.

c
c  read list of files and frequencies to be analyzed and file to output results
c  Usually will pipe in data from some file like fasearrayinp1


c#################################################################################
c######################### read the eqlistper file ################################
      read(*,*) nevents
      nobs = 0
      do iev = 1, nevents
        read(*,*) nsta(iev),idnum(iev)
        nobs = nobs+ 2*nsta(iev)
      do i = 1, nsta(iev)
      read(*, '(a)') fn(iev,i)
c     write(*, '(a)') fn(iev,i)
      enddo
      enddo
      write(*,*) 'it is ok there'
      read(*,*) nfreq
c always = 1
      read(*,*) (freq(j), j= 1, nfreq)
c frequency of interest
      read(*,*) foutput0
c detail.???.inp
      read(*,*) fsummary0
c summar.???.inp
      read(*,*) finvrsnodes
c gridnodes
      read(*,*) fftinput
c phampcor.???.inp
      read(*,*) fvariance0
c covar.???.inp
      read(*,*) dummy
c mavamp.???.inp
c      read(*,*) fresmax
      read(*,*) fvelarea0
c temp.???.inp
      read(*,*) dummy
c velarea.???.inp
      read(*,*) iterlimit, wlambda, dampvel,dampaniso
c four parameters
      read(*,*) unifvel
      read(*,'(a)') sensfn    
c sensitivity kernels

      iblank1 = blank(foutput0)
      iblank2 = blank(fsummary0)
      iblank3 = blank(fvariance0)
      iblank4 = blank(fvelarea0)

      foutput = foutput0(1:iblank1) // '.sa360kern'
      fsummary = fsummary0(1:iblank2) // '.sa360kern'
      fvariance = fvariance0(1:iblank3) //'.sa360kern'
      fvelarea = fvelarea0(1:iblank4) //'.sa360kern'
c
      open(10, file = foutput, status='new')
      open(11, file = fsummary,status='new')
      open(12, file = fftinput)
      open(14, file = "followit12")
      open(15, file = finvrsnodes)
      open(16, file = fvariance)
      open(66,file = sensfn)
c      open(80,file  = fresmax)



ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc find the station in the stationid.dat CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc find the station in the stationid.dat CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc this part need to be changed in order get station name 
 
      do ista = 1, maxnsta
        nevntsta(ista) = 0
      enddo

      do iev = 1, nevents
        do ista = 1, nsta(iev)
c          idtemp = index(fn(iev,ista),'..LHZ.00')
c          idstaname = fn(iev,ista)((idtemp-4):(idtemp-1))
          idstaname = fn(iev,ista)
          istanum(iev,ista) = findid(idstaname)
          nevntsta(istanum(iev,ista)) = nevntsta(istanum(iev,ista)) +1
        enddo
      enddo

c  count number of stations that have events and assign each station number a
c  variable number

cc  do following step to find the instrument type for each station
        
      do ista = 1, maxnsta
         istatype(ista) = findstatype(ista)
c        write(*,*) 'type', ista, istatype(ista)
      enddo

         do itype =1, maxnstatype
           itypeflag(itype) = 0
         enddo
         do i =1, maxnsta
           istaflag(i) = 0
         enddo

       do iev = 1, nevents
        do ista = 1, nsta(iev)
         
         do ijk = 1,maxnsta                     
            if(istanum(iev,ista) .eq. ijk ) then
              istaflag(ijk) = 1 
            endif              
         enddo
          do itype =1, maxnstatype
           if(istatype(istanum(iev,ista)) .eq. itype ) then
             itypeflag(itype) = 1 
           endif
          enddo
        enddo
       enddo

         ntype =0
         do itype =1, maxnstatype
            if( itypeflag(itype) .eq. 1) then
               ntype = ntype + 1
               iptype(itype) = ntype
               sttypename(ntype) = itype
            endif
         enddo

cc    nsta_real is the number of stations having data.          
         nsta_real = 0
         do ista =1, maxnsta
            if( istaflag(ista) .eq. 1) then
               nsta_real = nsta_real + 1
            endif
         enddo


       do iev = 1, nevents
        do ista = 1, nsta(iev)
          idtemp = index(fn(iev,ista),'..LHZ.00')
          idstaname = fn(iev,ista)((idtemp-4):(idtemp-1))

          if( istatype(istanum(iev,ista)) == 0 
     1   .or. iptype(istatype(istanum(iev,ista)))==0 ) then
         write(*,*) 'station type has problems'
           stop
        endif

c          write(*,*) iev,ista,istatype(istanum(iev,ista)),
c     1        iptype(istatype(istanum(iev,ista))),idstaname

        enddo
       enddo
       write(*,*) 'ntype',ntype,nsta_real

ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc CCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC CCCCCCCCCCCCCCCCCC




c################################################################
cc    read sensitivity kernels, the node file and data of phase and amplitude
c################################################################
          do ixkern = 1,maxnxints
        do iykern = 1,maxnyints
            sensitivity(ixkern,iykern) = 0.
                ampsens(ixkern,iykern) = 0.
        enddo
        enddo
      read(66,*) nxkern, xbegkern,dxkern
      read(66,*) nykern, ybegkern,dykern
          do ixkern = 1,nxkern
        do iykern = 1,nykern
         read(66,*) xdummy, ydummy,sensitivity(ixkern,iykern)
     1                                  ,ampsens(ixkern,iykern)
        enddo
        enddo
       
          write(*,*) 'sensitivity kernel'

c  read velocity nodes, initial velocities, box for limits of tomography
c  position in decimal degrees, phase velocity and aniso coeff in km/s, 
      read(15,*) dummy
      read(15,*) nnodes
      do i = 1,nnodes   
        read(15,*) nodelat(i),nodelon(i)
      enddo
      do i = 1, 4
        read(15, *) boxlat(i),boxlon(i)
      enddo
      read(15, *) ncol
      read(15,*) dxnode,dynode



      write(11,*) foutput

c  read data
      do iev = 1, nevents
        read(12,*) iv
        do ista = 1, nsta(iev)
          read(12,*) beg(iev,ista)
          read(12,*) stadist(iev,ista),staazi(iev,ista), bazi(iev,ista),
     1    stadelt(iev,ista), stalat(iev,ista), stalon(iev,ista)
          read(12,*) staamp(iev,ista,1),staph(iev,ista,1)

 
        enddo
      enddo

c################################################################
c################################################################




c  begin inversion loop over frequencies

      do ifreq = 1, nfreq

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      initialize the inverison by setting up parametres, coordinates, and so on.
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c  assume a priori data covariance = 0.1  Function of this constant value is
c  to make damping based on real estimates about the right size.

        do iev = 1, nevents
          stddevdata(iev) = 0.2
        enddo
        npnoamp = 6*nevents + 1*nnodes
        np = npnoamp + ntype*2
        kj = nnodes/ncol
        i6 = 6*nevents

c   ***************************************
c
c  annlscl is crudely equivalent to standard deviation for simulated annealing 
c  search for source parameters
        annlscl(1) = .25
        annlscl(2) = .25
        annlscl(3) = 15.*convdeg
        annlscl(4) = 15.*convdeg
        annlscl(5) = .2
        annlscl(6) = .2
        do iev = 1, nevents
          ip = (iev-1)*6
          covinv(1+ip) = 1./(0.40**2)
          covinv(2+ip) = 1./(0.40**2)
          covinv(3+ip) = 1./((10.*convdeg)**2)
          covinv(4+ip) = 1./((10.*convdeg)**2)
          covinv(5+ip) = 1./(.25**2)
          covinv(6+ip) = 1./(.25**2)

          startamp1(iev) = 0.5
          startamp2(iev) = 0.5
          stazim1(iev) = -7.0*convdeg
          stazim2(iev) = 7.0*convdeg
          stphase1(iev) = 0.0
          stphase2(iev) = 0.0
          origmod(1+ip) = startamp1(iev)
          origmod(2+ip) = startamp2(iev)
          origmod(3+ip) = stazim1(iev)
          origmod(4+ip) = stazim2(iev)
          origmod(5+ip) = stphase1(iev)
          origmod(6+ip) = stphase2(iev)
          crrntmod(1+ip) = startamp1(iev)
          crrntmod(2+ip) = startamp2(iev)
          crrntmod(3+ip) = stazim1(iev)
          crrntmod(4+ip) = stazim2(iev)
          crrntmod(5+ip) = stphase1(iev)
          crrntmod(6+ip) = stphase2(iev)
        enddo
        do ii= 1, nnodes
          ip = i6 + ii
          origmod(ip) = unifvel
          covinv(ip) = 1./(dampvel**2)
          crrntmod(ip) = origmod(ip)
        enddo

        do ii = 1, ntype
          ampmult(ii) = 1.0
          ip = npnoamp+ii
          origmod(ip) = ampmult(ii)
          crrntmod(ip) = origmod(ip)
          covinv(ip) = 1.0/(0.30**2)
        enddo

        do ii = 1, ntype
          phcor(ii) = 0.
          ip = npnoamp + ntype + ii
          origmod(ip) = phcor(ii)
          crrntmod(ip) = origmod(ip)
          covinv(ip) = 1.0/(0.1**2)
        enddo


cc^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c  increase the variance for edges for velocity 
cc^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        varfac2 = 10.
        do itp = 1,1
         ityp = nnodes*(itp-1)
c  right end 
         do ijk = i6+1,i6+kj
           jk = ijk+ityp
           covinv(jk) = covinv(jk)/varfac2
         enddo
c  left end
         do ijk = i6+nnodes-kj+1,i6+nnodes
           jk = ijk+ityp
           covinv(jk) = covinv(jk)/varfac2
         enddo
c  top
         do ijk = i6+kj+1,i6+nnodes-2*kj+1,kj
           jk = ijk+ityp
           covinv(jk) = covinv(jk)/varfac2
         enddo
c  bottom
         do ijk = i6+2*kj,i6+nnodes-kj,kj
           jk = ijk+ityp
           covinv(jk) = covinv(jk)/varfac2
         enddo
        enddo
cc^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
cc^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      do iev = 1, nevents
c
c  Also find normalizing amplitude to equalize earthquakes of different
c  size, using rms amplitude (old version used largest amplitude)

        diffamp = 10000000000.
        amprms(iev,ifreq) = 0.0
        iref(iev) = 1
      stasmall = 100000.
        do ista = 1, nsta(iev)
              amprms(iev,ifreq) = amprms(iev,ifreq)
     1                            + staamp(iev,ista,ifreq)**2
        enddo

        amprms(iev,ifreq) = sqrt(amprms(iev,ifreq)/nsta(iev))

        do ista = 1, nsta(iev)
           if (abs(staamp(iev,ista,ifreq)-amprms(iev,ifreq))
     1             .lt. diffamp) then
            diffamp = abs(staamp(iev,ista,ifreq)-amprms(iev,ifreq))
            iref(iev) = ista
           endif
        enddo
c  
c  Use reference station for to set up local coordinate system
c  based on reference station at zero, zero.  Use distances to stations rather
c  than absolute coordinates as way of correcting for curving wavefront.  This
c  will tend to favor keeping azimuth at original azimuth along great circle.
c  + x direction is in direction of propagation along great circle path.
c  +y is 90 deg counterclockwise from x
c   staazi are measured clockwise from north


      xsta(iev,iref(iev))  = 0.0
      ysta(iev,iref(iev))  = 0.0
      rloc(iev,iref(iev))  = 0.0
      azloc(iev,iref(iev)) = 0.0

      do ista = 1, nsta(iev)
        if (ista.ne.iref(iev)) then
        xsta(iev,ista) = stadist(iev,ista) - stadist(iev,iref(iev))
        azidiff = staazi(iev,iref(iev)) - staazi(iev,ista)
        if (azidiff.gt.180.) azidiff = azidiff -360.
        if (azidiff.lt.-180.) azidiff = azidiff +360.
        ysta(iev,ista) = circ*sin(stadelt(iev,ista)*convdeg)*
     1                   azidiff
        rloc(iev,ista) = sqrt(xsta(iev,ista)*xsta(iev,ista) + 
     1                        ysta(iev,ista)*ysta(iev,ista))
        azloc(iev,ista) = atan2(ysta(iev,ista),xsta(iev,ista))
        endif
      enddo

c  calculate apparent pole position of earthquake for spherical earth
c  instead of real location so that coordinate system for nodes will agree
c  with that for stations 
      call gohead(stalat(iev,iref(iev)),stalon(iev,iref(iev)),
     1    stadelt(iev,iref(iev)),bazi(iev,iref(iev)),applat,applon)
      call disthead(applat,applon,stalat(iev,iref(iev)),
     1    stalon(iev,iref(iev)),delta,tazimref)
c  now calculate x,y of each node 
c  this approach seems to work pretty well - agrees within about 0.1% with
c  station calculations based on elliptical great circle distance in terms
c  of relative position
      appcirc = stadist(iev,iref(iev))/stadelt(iev,iref(iev))
      do inode = 1, nnodes
        call disthead(applat,applon,nodelat(inode),nodelon(inode)
     1                                       ,delta,tazim)
        call disthead(nodelat(inode),nodelon(inode),applat,applon
     1                                       ,delta,bazimnd(inode))
        xnode(iev,inode) = appcirc*delta - stadist(iev,iref(iev))
        azidiff = tazimref - tazim
        if (azidiff.gt.180.) azidiff = azidiff -360.
        if (azidiff.lt.-180.) azidiff = azidiff +360.
        ynode(iev,inode) = appcirc*sin(delta*convdeg)*azidiff
      enddo
c  similarly for outlines of region of interest
      do ibox = 1, 4
        call disthead(applat,applon,boxlat(ibox),boxlon(ibox)
     1                                       ,delta,tazim)
        xbox(iev,ibox) = appcirc*delta - stadist(iev,iref(iev))
        azidiff = tazimref - tazim
        if (azidiff.gt.180.) azidiff = azidiff -360.
        if (azidiff.lt.-180.) azidiff = azidiff +360.
        ybox(iev,ibox) = appcirc*sin(delta*convdeg)*azidiff
c      write(*,*) ibox, xbox(ibox),ybox(ibox),delta, tazim
      enddo

              
c  generate real and imaginary components normalized by amplitude at reference
c  station and compared to phase at reference station
c  Phase shift relative to reference corrected for any difference in start time
                  
        do ista = 1, nsta(iev)

          dph(iev,ista) = staph(iev,ista,ifreq)-staph(iev,iref(iev),
     1      ifreq) + freq(ifreq)*(beg(iev,ista)-beg(iev,iref(iev)))
          streal(iev,ista,ifreq)=staamp(iev,ista,ifreq)
     1            *cos(dph(iev,ista)*twopi)/amprms(iev,ifreq)
          stimag(iev,ista,ifreq)=-staamp(iev,ista,ifreq)
     1            *sin(dph(iev,ista)*twopi)/amprms(iev,ifreq)

        enddo           
      enddo       

 
cc  below is to find the mininum point of the four box points
        do 123 iev = 1, nevents
        do 124 ideg = 1,ndeg
         tpstazim = ((ideg-1.) - (ndeg-1)/2)*convdeg

c  find minimum xbox
       iflag = 1 
       templen0 =   
     1   xbox(iev,1)*cos(tpstazim)
     1 + ybox(iev,1)*sin(tpstazim) 

      do ibox = 1,4

      xmintemp = 
     1   xbox(iev,ibox)*cos(tpstazim)
     1 + ybox(iev,ibox)*sin(tpstazim) 

         if ( xmintemp .lt. templen0) then
            templen0 = xmintemp
            iflag = ibox
       endif 
      enddo

      xmin(iev,ideg) = 
     1   xbox(iev,iflag)*cos(tpstazim)
     1 + ybox(iev,iflag)*sin(tpstazim) 

  124    enddo
  123    enddo

       write(*,*)  'before inversion '
         write(*,*) 'nobs' ,nobs, 'np ',np



c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


c
c  *************************
c  iterate from here
c  *************************
        iter = 1
        icnt = 1
      nobs = nobs + 2 


100     continue

       do iev = 20, 20
        do ista = 1, nsta(iev)

c          write(*,*) iev,ista,istatype(istanum(iev,ista)),
c     1        iptype(istatype(istanum(iev,ista)))

        enddo
       enddo

        dstacnt = ntype
        dstaph  = 0. 

      write(*,*) iter,icnt,ntype


cc@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c            g for station correction terms
cc@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        sumampcor = 0.0
        do ii = 1, ntype
          sumampcor = sumampcor + ampmult(ii)
        enddo
        d(nobs-1) = (dstacnt - sumampcor)/0.2        


        sumphcor = 0.0
        do ii = 1, ntype
          sumphcor = sumphcor + phcor(ii)
        enddo
        d(nobs) = (dstaph - sumphcor)/0.2        

      write(*,*) 'before g'

         allocate(g(nobs,np))

        do irow = 1, nobs
        do icol = 1, np
          g(irow,icol) = 0.0
        enddo
        enddo

c xx
        do ii = 1,ntype
          icol  = npnoamp+ii
          icol2 = npnoamp+ii + ntype
          g(nobs-1,icol) =  1.0/0.2
          g(nobs  ,icol2) = 1.0/0.2
        enddo 

cc@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cc@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



c  begin loop over events for residuals and partial derivatives

        naddat = 0
        write(*,*) 'before setup wgt'

        allocate(wgtnode1(nsta_real,nnodes,ndeg))
        allocate(ampwgtnode1(nsta_real,nnodes,ndeg) ) 

        call system('date')

cc    Go through each event

        do 60  iev = 1, nevents

c  calculate current apparent velocity at each node
        do ii = 1, nnodes
          iii = ii+i6
          jjj = i6+nnodes+(ii-1)/kj + 1
          jjjj = jjj + iages
          appvel(ii) = crrntmod(iii)
        enddo


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
cc   calculate the sensitivity kernel and the contribution to amplitude and phase from sensitivity 
cc      for each possible angle bwtween -40 to 40
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
        do 125 ideg = 1,ndeg
         tpstazim = ((ideg-1.) - (ndeg-1)/2)*convdeg
        do 102 ista = 1, nsta(iev)

        do ii = 1, nnodes
            wgtnode1(ista,ii,ideg) = 0.0
         ampwgtnode1(ista,ii,ideg) = 0.0
        enddo

       xstatemp = 
     1   xsta(iev,ista)*cos(tpstazim)
     1 + ysta(iev,ista)*sin(tpstazim)               
     
       ystatemp = 
     1  -xsta(iev,ista)*sin(tpstazim)
     1 + ysta(iev,ista)*cos(tpstazim)     

   
c  normalize weights to be equivalent to distance

         senssum1 = 0. 

            do ii = 1,nnodes

       xnodetemp = 
     1   xnode(iev,ii)*cos(tpstazim)
     1 + ynode(iev,ii)*sin(tpstazim)               
     
       ynodetemp = 
     1  -xnode(iev,ii)*sin(tpstazim)
     1 + ynode(iev,ii)*cos(tpstazim)     


         xstanode = xnodetemp - xstatemp
       ystanode = ynodetemp - ystatemp

       if( xnodetemp .ge. xmin(iev,ideg)) then
          ixindex = int( xstanode/dxkern ) + (nxkern+1)/2
        iyindex = int( ystanode/dykern ) + (nykern+1)/2

          if(  ixindex .lt.1 .or. ixindex .gt. nxkern
     1    .or. iyindex .lt.1 .or. iyindex .gt. nxkern) then 


c         write(*,*) 'ixindex,iyindex', ixindex,iyindex
           wgtnode1(ista,ii,ideg)    = 0.
           ampwgtnode1(ista,ii,ideg) = 0.  

         else

              wgtnode1(ista,ii,ideg) = sensitivity(ixindex,iyindex) 
     1                    *(dxnode*dynode)/(dxkern*dykern)


           ampwgtnode1(ista,ii,ideg) =     ampsens(ixindex,iyindex)
     1                    *(dxnode*dynode)/(dxkern*dykern)
           endif 

      else 
           wgtnode1(ista,ii,ideg)    = 0.
           ampwgtnode1(ista,ii,ideg) = 0.  
        endif

          enddo


       dphase(ista,ideg) = 0.
        dampper(ista,ideg) = 0.

          do inode = 1, nnodes
          dphase(ista,ideg) = dphase(ista,ideg) 
     1  + (1.0/twopi)*wgtnode1(ista,inode,ideg)
     1           *(appvel(inode)-unifvel)/unifvel

          enddo


        do inode =1, nnodes
          dampper(ista,ideg) = dampper(ista,ideg) 
     1  +ampwgtnode1(ista,inode,ideg)*
     1   (appvel(inode)-unifvel)/unifvel

          enddo

 102    enddo
 125     enddo

cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

cc end of kernels



ccc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ccc +++++++++++++++++++ begin annealing ++++++++++++++++++++++++++++
         ip = (iev-1)*6

c  start simulated annealing setup
        nrestart = 3
        bestfit = 1.0E+16
        irestart = 0
c  set up original simplex of points to start simulated annealing process
c  starting with current model and one point perturbed in direction of 
c  each variable
        do ismplx = 1, 7
           iprtrb = ismplx-1
           do ismp = 1,6
             p(ismplx,ismp) = crrntmod(ip+ismp)
             if (ismp.eq.iprtrb) p(ismplx,ismp) = crrntmod(ip+ismp) +
     1                                             annlscl(ismp)
             pv(ismp) = p(ismplx,ismp)
c          write(*,*) pv(ismp)
           enddo
c  calculate misfit at each vertex
           smsft(ismplx) = misfit(pv)
c        write(*,*) misfit(pv)
        enddo

c  set starting temperature, tolerance, cooling factor, loops, and iteration counter
c        write(14,*)
c        write(14,*) 'event',iev

498   temptr = 100.
        tempinit = temptr
        cfac = 0.5
        nloop = 14
        niter = 100
        ftol = .0001
        loopann = 0
499     iterann = niter
        call anneal(p,smsft,7,6,6,pb,bestfit,ftol,iterann,temptr)
        loopann = loopann + 1
        if (loopann.le.nloop) then
           temptr = temptr*cfac
           go to 499
        endif
        irestart= irestart +1
c  use three other starting models to make sure that best model is not missed
c  due to difficulty of finding absolute minimum
          if (irestart.le.nrestart) then
            if ((irestart.eq.1).or.(irestart.eq.2)) then
              pppb(1) = 0.7
              pppb(2) = 0.3
              pppb(5) = 0.0
              pppb(6) = 0.0      
              pppb(3) = 0.*convdeg
              pppb(4) = -10.*convdeg
            endif
            if (irestart.eq.1) pppb(4) = 10.*convdeg
            if (irestart.eq.3) then
              pppb(1) = 0.5
              pppb(2) = 0.5
              pppb(5) = 0.0
              pppb(6) = 0.0      
              pppb(3) = 7.*convdeg
              pppb(4) = -7.*convdeg
            endif
            do ismplx = 1, 7
              iprtrb = ismplx-1
              do ismp = 1,6
                p(ismplx,ismp) = pppb(ismp)
                if (ismp.eq.iprtrb) p(ismplx,ismp) = pppb(ismp) +
     1                                             annlscl(ismp)
                pv(ismp) = p(ismplx,ismp)
              enddo
c  calculate misfit at each vertex
              smsft(ismplx) = misfit(pv)
c              write(*,*)      misfit(pv)
            enddo
           
            go to 498
         endif 

c           write(*,*) misfit(pb)

c  use best model as event starting model for linearized inversion
        startamp1(iev) = pb(1)
        startamp2(iev) = pb(2)
        stazim1(iev) = pb(3)
        stazim2(iev) = pb(4)
        stphase1(iev) = pb(5)
        stphase2(iev) = pb(6)
cc+++++++++++++++++++++++ end of annealing++++++++++++++++++++++++++++++++++
ccc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


         write(*,*) 'after annealing',iev, 
     1     sqrt(misfit(pb)/nsta(iev))*stddevdata(iev)
c
      write(*,*) startamp1(iev),stazim1(iev)/convdeg,stphase1(iev)
      write(*,*) startamp2(iev),stazim2(iev)/convdeg,stphase2(iev) 


       if(stazim1(iev) .lt. -38*convdeg) stazim1(iev) = -38*convdeg
       if(stazim1(iev) .gt. 38*convdeg)  stazim1(iev) = 38*convdeg
       if(stazim2(iev) .lt. -38*convdeg) stazim2(iev) = -38*convdeg
       if(stazim2(iev) .gt. 38*convdeg)  stazim2(iev) = 38*convdeg


         ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2) + 1
         ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2) + 1



c  find minimum xbox of the grid nodes for each event

       do 10 ista = 1, nsta(iev)

       xstatemp = 
     1   xsta(iev,ista)*cos(stazim1(iev))
     1 + ysta(iev,ista)*sin(stazim1(iev))               
     
       ystatemp = 
     1  -xsta(iev,ista)*sin(stazim1(iev))
     1 + ysta(iev,ista)*cos(stazim1(iev))     

      dtime1(ista) = (xstatemp-xmin(iev,ideg1))*(1.0/unifvel)

  10    enddo


c       calculate for second plane wave

        do 20 ista = 1, nsta(iev)
         
       xstatemp = 
     1   xsta(iev,ista)*cos(stazim2(iev))
     1 + ysta(iev,ista)*sin(stazim2(iev))               
     
       ystatemp = 
     1  -xsta(iev,ista)*sin(stazim2(iev))
     1 + ysta(iev,ista)*cos(stazim2(iev))     

      dtime2(ista) = (xstatemp-xmin(iev,ideg2))*(1.0/unifvel)

  20    enddo

      


c  calculate apparent average slowness or delay from edge of study area to 
c  each station

        do ista = 1, nsta(iev)
         phase1(iev,ista) = dtime1(ista)*freq(1)
         phase2(iev,ista) = dtime2(ista)*freq(1)
        dphase1(iev,ista) =   dphase(ista,ideg1)
        dphase2(iev,ista) =   dphase(ista,ideg2)  
          damp1per(iev,ista) = dampper(ista,ideg1)
        damp2per(iev,ista) = dampper(ista,ideg2) 
       enddo



c###########################################################################################################
c#calulate the change of ph1,ph2, damp1,damp2 with respect to azimuth at the reference station######## ista
c###########################################################################################################


      aziminc = 2.*convdeg
      stazim1(iev) = stazim1(iev) + aziminc
      stazim2(iev) = stazim2(iev) + aziminc


       if(stazim1(iev) .lt. -38*convdeg) stazim1(iev) = -38*convdeg
       if(stazim1(iev) .gt. 38*convdeg)  stazim1(iev) = 38*convdeg
       if(stazim2(iev) .lt. -38*convdeg) stazim2(iev) = -38*convdeg
       if(stazim2(iev) .gt. 38*convdeg)  stazim2(iev) = 38*convdeg

 
         ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2) + 1
         ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2) + 1



c  initialize weights

c      calculate for first plane wave


       xstatemp = 
     1   xsta(iev,iref(iev))*cos(stazim1(iev))
     1 + ysta(iev,iref(iev))*sin(stazim1(iev))               
     
       ystatemp = 
     1  -xsta(iev,iref(iev))*sin(stazim1(iev))
     1 + ysta(iev,iref(iev))*cos(stazim1(iev))     

   
c  normalize weights to be equivalent to distance


       dtime1ref = (xstatemp-xmin(iev,ideg1))*(1.0/unifvel)


c       calculate for sencond plane wave ist
c  initialize weights

       xstatemp = 
     1   xsta(iev,iref(iev))*cos(stazim2(iev))
     1 + ysta(iev,iref(iev))*sin(stazim2(iev))               
     
       ystatemp = 
     1  -xsta(iev,iref(iev))*sin(stazim2(iev))
     1 + ysta(iev,iref(iev))*cos(stazim2(iev))     

   
c  normalize weights to be equivalent to distance

            
        dtime2ref = (xstatemp-xmin(iev,ideg2))*(1.0/unifvel)

          

          phase1ref = dtime1ref*freq(1)
        phase2ref = dtime2ref*freq(1)

          dphase1ref =   dphase(iref(iev),ideg1)
        dphase2ref =   dphase(iref(iev),ideg2) 
          damp1perref = dampper(iref(iev),ideg1)
        damp2perref = dampper(iref(iev),ideg2) 

      stazim1(iev) = stazim1(iev) - aziminc
      stazim2(iev) = stazim2(iev) - aziminc

c################### end of for refrence station #################




        do 50 ista = 1,nsta(iev)

           iptptemp = iptype(istatype(istanum(iev,ista)))
           ampadj = ampmult(iptptemp)
           phadj =   phcor(iptptemp)

        prefase1 = ( (phase1(iev,ista)+dphase1(iev,ista))
     1  - (phase1(iev,iref(iev))+dphase1(iev,iref(iev))))
     1             + stphase1(iev)
        prefase2 = ( (phase2(iev,ista)+dphase2(iev,ista))
     1  - (phase2(iev,iref(iev))+dphase2(iev,iref(iev))))
     1             + stphase2(iev)

          prefase1 = prefase1 + phadj
          prefase2 = prefase2 + phadj

        staamp1 = startamp1(iev)*(1.+damp1per(iev,ista))
        staamp2 = startamp2(iev)*(1.+damp2per(iev,ista))


          cosph1 = cos(prefase1*twopi)
          cosph2 = cos(prefase2*twopi)
          sinph1 = sin(prefase1*twopi)
          sinph2 = sin(prefase2*twopi)
          
          prereal = staamp1*cosph1+staamp2*cosph2  
        preimag = -1.0*(staamp1*sinph1+ staamp2*sinph2)

          prereal = prereal*ampadj
          preimag = preimag*ampadj


          kreal = ista + naddat
          kimag = kreal + nsta(iev)
          d(kreal) = (streal(iev,ista,ifreq) - prereal)/stddevdata(iev)
          d(kimag) = (stimag(iev,ista,ifreq) - preimag)/stddevdata(iev)


c###########################################################################################################
c                 calculate parph1azim, parph2azim, pardamp1azim,pardamp2azim
c############################ calculate parph1azim,parph2azim, #############################################


c  number of intervals = nints, interval length = actxint
      aziminc = 2.*convdeg  
      stazim1(iev) = stazim1(iev) + aziminc
      stazim2(iev) = stazim2(iev) + aziminc


       if(stazim1(iev) .lt. -38*convdeg) stazim1(iev) = -38*convdeg
       if(stazim1(iev) .gt.  38*convdeg) stazim1(iev) = 38*convdeg
       if(stazim2(iev) .lt. -38*convdeg) stazim2(iev) = -38*convdeg
       if(stazim2(iev) .gt.  38*convdeg) stazim2(iev) = 38*convdeg

         ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2) + 1
         ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2) + 1

c  initialize weights

c      calculate for first plane wave


       xstatemp = 
     1   xsta(iev,ista)*cos(stazim1(iev))
     1 + ysta(iev,ista)*sin(stazim1(iev))               
     
       ystatemp = 
     1  -xsta(iev,ista)*sin(stazim1(iev))
     1 + ysta(iev,ista)*cos(stazim1(iev))     

   
c  normalize weights to be equivalent to distance




         dtime1temp = (xstatemp-xmin(iev,ideg1))*(1.0/unifvel)
 


c       calculate for sencond plane wave
         
c  initialize weights


       xstatemp = 
     1   xsta(iev,ista)*cos(stazim2(iev))
     1 + ysta(iev,ista)*sin(stazim2(iev))               
     
       ystatemp = 
     1  -xsta(iev,ista)*sin(stazim2(iev))
     1 + ysta(iev,ista)*cos(stazim2(iev))     

   
c  normalize weights to be equivalent to distance

            

        dtime2temp = (xstatemp-xmin(iev,ideg2))*(1.0/unifvel)

          phase1temp = dtime1temp*freq(1)
        phase2temp = dtime2temp*freq(1)
        dphase1temp =   dphase(ista,ideg1)
        dphase2temp =   dphase(ista,ideg2)
          damp1pertemp = dampper(ista,ideg1)
        damp2pertemp = dampper(ista,ideg2) 
  
      stazim1(iev) = stazim1(iev) - aziminc
      stazim2(iev) = stazim2(iev) - aziminc


       if(stazim1(iev) .lt. -38*convdeg) stazim1(iev) = -38*convdeg
       if(stazim1(iev) .gt.  38*convdeg) stazim1(iev) = 38*convdeg
       if(stazim2(iev) .lt. -38*convdeg) stazim2(iev) = -38*convdeg
       if(stazim2(iev) .gt.  38*convdeg) stazim2(iev) = 38*convdeg

         ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2) + 1
         ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2) + 1


c  normalize weights to be equivalent to distance
          parph1azim = ( (phase1temp+dphase1temp)-
     1          (phase1(iev,ista)+dphase1(iev,ista)))/aziminc
     1                 - ( (phase1ref+dphase1ref) -
     1          (phase1(iev,iref(iev))+dphase1(iev,iref(iev))))/aziminc

          parph2azim = ( (phase2temp+dphase2temp)-
     1          (phase2(iev,ista)+dphase2(iev,ista)))/aziminc
     1                 - ( (phase2ref+dphase2ref) -
     1          (phase2(iev,iref(iev))+dphase2(iev,iref(iev))))/aziminc

          pardamp1azim = (damp1pertemp - damp1per(iev,ista))/aziminc
        pardamp2azim = (damp2pertemp - damp2per(iev,ista))/aziminc

c########## end of parph1azim,parph2azim,pardamp1azim,pardamp2azim######## 
c###########################################################################################################
c###########################################################################################################


c  partial derivatives for station amplitude correction factors and attenuation
          g(kreal,npnoamp+iptptemp) = 
     1      prereal/ampadj/stddevdata(iev)
          g(kimag,npnoamp+iptptemp) = 
     1      preimag/ampadj/stddevdata(iev)

c  partial derivatives for station phase correction

          g(kreal,npnoamp+ntype+iptptemp) =      
     1             (-startamp1(iev)*sinph1*twopi
     2           -startamp2(iev)*sinph2*twopi)*ampadj
     3                  /stddevdata(iev)
                         
          g(kimag,npnoamp+ntype+iptptemp) =
     1             (-startamp1(iev)*cosph1*twopi
     2           -startamp2(iev)*cosph2*twopi)*ampadj
     3                   /stddevdata(iev)



       do ii = 1, nnodes

       parph1v = (1.0/twopi)*wgtnode1(ista,ii,ideg1)/unifvel
     1         - (1.0/twopi)*wgtnode1(iref(iev),ii,ideg1)/unifvel

       parph2v =  (1.0/twopi)*wgtnode1(ista,ii,ideg2)/unifvel
     1          - (1.0/twopi)*wgtnode1(iref(iev),ii,ideg2)/unifvel

       paramp1v = startamp1(iev)*ampwgtnode1(ista,ii,ideg1)/unifvel 
       paramp2v = startamp2(iev)*ampwgtnode1(ista,ii,ideg2)/unifvel

c  partial derivatives with respect to velocity, & cos2theta &sin2theta
c***************************** changed for no aniso
c            jjjage = (ii-1)/kj +1 + i6 + nnodes
c            jjjjage = jjjage + iages

        g(kreal,i6+ii)= (-startamp1(iev)*(1.0+damp1per(iev,ista))*
     1                    twopi*sinph1*parph1v
     1                   -startamp2(iev)*(1.0+damp2per(iev,ista))*
     1                     twopi*sinph2*parph2v 
     3                 + paramp1v*cosph1 + paramp2v*cosph2 ) 
     2                *ampadj/stddevdata(iev) 

        g(kimag,i6+ii)= (-startamp1(iev)*(1.0+damp1per(iev,ista))*
     1                  twopi*cosph1*parph1v
     1                 -startamp2(iev)*(1.0+damp2per(iev,ista))*
     1                          twopi*cosph2*parph2v
     3              - (paramp1v*sinph1  +  paramp2v*sinph2 ) )
     2                *ampadj/stddevdata(iev) 

       enddo
c  partial derivatives in order are with respect to amplitudes,
c  azimuths, starting phases and  slowness
          ip = (iev-1)*6
        g(kreal,1+ip) = (1.0+damp1per(iev,ista))*cosph1
     1                     *ampadj/stddevdata(iev)
        g(kreal,2+ip) = (1.0+damp2per(iev,ista))*cosph2
     1                     *ampadj/stddevdata(iev)

        g(kreal,3+ip) = (-startamp1(iev)*(1.0+damp1per(iev,ista))*sinph1
     1                      *parph1azim*twopi
     1                    + startamp1(iev)*pardamp1azim*cosph1)
     1                      *ampadj/stddevdata(iev)

        g(kreal,4+ip) = (-startamp2(iev)*(1.0+damp2per(iev,ista))*sinph2
     1                      *parph2azim*twopi
     1                    + startamp2(iev)*pardamp2azim*cosph2) 
     1                     *ampadj/stddevdata(iev)

        g(kreal,5+ip) = -startamp1(iev)*(1.0+damp1per(iev,ista))*sinph1
     1                      *twopi
     1                     *ampadj/stddevdata(iev)

        g(kreal,6+ip) = -startamp2(iev)*(1.0+damp2per(iev,ista))*sinph2
     1                      *twopi
     1                      *ampadj/stddevdata(iev)

        g(kimag,1+ip) = -(1.0+damp1per(iev,ista))*sinph1
     1                      *ampadj/stddevdata(iev)

        g(kimag,2+ip) = -(1.0+damp2per(iev,ista))*sinph2
     1                      *ampadj/stddevdata(iev)

        g(kimag,3+ip) = -(startamp1(iev)*(1.0+damp1per(iev,ista))*cosph1
     1                      *parph1azim*twopi
     1                      + startamp1(iev)*pardamp1azim*sinph1)
     1                      *ampadj/stddevdata(iev)

        g(kimag,4+ip) = -(startamp2(iev)*(1.0+damp2per(iev,ista))*cosph2
     1                      *parph2azim*twopi
     1                      + startamp2(iev)*pardamp2azim*sinph2)
     1                      *ampadj/stddevdata(iev)

        g(kimag,5+ip) = -startamp1(iev)*(1.0+damp1per(iev,ista))*cosph1
     1                      *twopi
     1                     *ampadj/stddevdata(iev)

        g(kimag,6+ip) = -startamp2(iev)*(1.0+damp2per(iev,ista))*cosph2
     1                      *twopi
     1                      *ampadj/stddevdata(iev)
c  normalize by stddevdata
 50    enddo

          naddat = naddat + 2*nsta(iev)
 60     enddo



        write(*,*) "end events for partial derivatives"
        call system('date')
         write(*,*) 'nobs' ,nobs, 'np ',np

        deallocate(wgtnode1,ampwgtnode1)


c
c  Calculate gtg and gtd
c


        allocate(gtg(np,np))

        do j = 1, np
c          write(*,*) j
          gtd(j) = 0.0
          do i = 1, nobs
            gtd(j) = gtd(j) + g(i,j)*d(i)
          enddo
c   add to gtd Tarantola term penalizing misfit to original starting model
c   but skip for wave parameters
c
          if (j.le.i6) gtdcmm(j) = gtd(j)       
          if (j.gt.i6) gtdcmm(j) = gtd(j) - 
     1                       covinv(j)*(crrntmod(j)-origmod(j))
c  construct gtg  
          do jj = 1,j
            gtg(jj,j) = 0.0
            do i = 1, nobs
              gtg(jj,j)= gtg(jj,j) + g(i,jj)*g(i,j)
            enddo
            gtg(j,jj) = gtg(jj,j)
cc            savegtg(j,jj) = gtg(j,jj)
cc            savegtg(jj,j) = gtg(jj,j)
          enddo
          gtg(j,j) = gtg(j,j) + covinv(j)
        enddo

        write(*,*) 'deallocate(g)'

         deallocate(g)

        write(*,*) 'before matrix inversion'
        write(*,*) 'np  ',np,nparam



c  Invert gtg.  gtg will be destroyed.  
C  Not the most efficient approach because doesn't take advantage of 
c  symmetry of gtg.  Use LU decomposition from Press et al.


        call system('date')



c  Find change to starting model

        do i= 1, np
            change(i) = gtdcmm(i)
        enddo


        write(*,*) 'before dludcmp'
        call dludcmp(gtg,np,np,indx,ddd)
        write(*,*) 'after dludcmp'     
        call dlubksb(gtg,np,np,indx,change)
        write(*,*) 'after matrix inversion'
        call system('date')



        if( (icnt .eq. iterlimit) .and. (iter .eq. 2) ) then
          allocate(gtginv(np,np))
          do i= 1,np
           do j = 1, np
             gtginv(i,j)= 0.0D0
           enddo
             gtginv(i,i) =1.0D0
          enddo
          do j = 1,np
           call dlubksb(gtg,np,np,indx,gtginv(1,j))
          enddo
        endif 

        deallocate(gtg)




c  Find rank (sum of diagonals of resolution matrix), i.e., number of
c  pieces of information or number of independent model parameters
c  rank1 is contribution from source wave terms, rank2 from velocity
c  variables



 


       



cUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
c  Update current model
cUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
 
        naddat = 0
        do iev = 1,nevents
          ip = (iev-1)*6
          startamp1(iev) = startamp1(iev) + change(1+ip)
          startamp2(iev) = startamp2(iev) + change(2+ip)
          stazim1(iev) = stazim1(iev) + change(3+ip)
          stazim2(iev) = stazim2(iev) + change(4+ip)

       if(stazim1(iev) .lt. -38*convdeg) stazim1(iev) = -38*convdeg
       if(stazim1(iev) .gt. 38*convdeg)  stazim1(iev) = 38*convdeg
       if(stazim2(iev) .lt. -38*convdeg) stazim2(iev) = -38*convdeg
       if(stazim2(iev) .gt. 38*convdeg)  stazim2(iev) = 38*convdeg

          stphase1(iev) = stphase1(iev) + change(5+ip)
          stphase2(iev) = stphase2(iev) + change(6+ip)
          crrntmod(1+ip) = startamp1(iev)
          crrntmod(2+ip) = startamp2(iev)
          crrntmod(3+ip) = stazim1(iev)
          crrntmod(4+ip) = stazim2(iev)
          crrntmod(5+ip) = stphase1(iev)
          crrntmod(6+ip) = stphase2(iev)

        sumsq = 0.0
        sumsqph = 0.0

          do ista = 1,nsta(iev)
            dresid1=d(ista+naddat)*stddevdata(iev)
            dresid2=d(ista+naddat+nsta(iev))*stddevdata(iev)
            predamp = amprms(iev,ifreq)
     1        *sqrt((streal(iev,ista,ifreq)-dresid1)**2 
     2                + (stimag(iev,ista,ifreq)-dresid2)**2)
            prefase = atan2(-(stimag(iev,ista,ifreq)-dresid2),
     1                (streal(iev,ista,ifreq)-dresid1))/twopi
           if (prefase - dph(iev,ista).gt.0.5) prefase = prefase-1.0
           if (prefase - dph(iev,ista).lt.-0.5) prefase = prefase+1.0   

           sumsq = sumsq+(d(ista+naddat)**2+d(ista+naddat+nsta(iev))**2)
     1               *stddevdata(iev)**2
           sumsqph = sumsqph + (dph(iev,ista)-prefase)**2
          enddo

             rmsph = sqrt(sumsqph/nsta(iev))/freq(ifreq)
             sigma2 = sqrt(sumsq/(2*nsta(iev)-8))



          naddat = naddat + 2*nsta(iev)

c     write(10,*) 'iev', idnum(iev),'  sigma2 ',sigma2,  'rmsph ',rmsph
c     write(10,*) stazim1(iev)/convdeg,startamp1(iev),stphase1(iev)
c     write(10,*) stazim2(iev)/convdeg,startamp2(iev),stphase2(iev)

        enddo




        do ii = 1, nnodes
          crrntmod(ii+i6) = crrntmod(ii+i6) + change(ii+i6)
        enddo


        do ii = 1, ntype
          ip = ii + npnoamp
          ampmult(ii) = ampmult(ii) + change(ip)
          crrntmod(ip) = ampmult(ii)
        enddo

        do ii = 1, ntype
          ip = ii + ntype + npnoamp
          phcor(ii) = phcor(ii) + change(ip)
          crrntmod(ip) = phcor(ii)
        enddo

        do ii = 1, ntype
       write(*,167)  sttypename(ii), ampmult(ii), phcor(ii) 
        enddo
cUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
cUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU


        chmax = 0.0
        do ii = 1,nnodes
          ip = i6 + ii
          chmax = dmax1(chmax, dabs(change(ip)))
        enddo
        write(14,*) icnt, chmax
        icnt = icnt +1
c  *******************************************
c  test for convergence and finish iteration
c  ********************************************
c        if ((chmax.gt.0.0001).and.(icnt.lt.500)) go to 100

        if ( icnt .eq. 3 .and. iter .eq. 1) goto 400

cc        if ((chmax.gt.0.0005).and ( icnt.le.iterlimit)) go to 100

         if ( ( icnt.le.iterlimit ) ) go to 100





c  The first time through, solutions are damped with data variance as an
c  estimate.  Go through again with a posteriori estimate of data variance
c  different for each event.
 400    continue
        if (iter.eq.1) then
           icnttot = icnt
           iter = 2
           minstd =.03
           naddat = 0
           do iev = 1, nevents
             sumsq = 0.0
             do ista = 1, 2*nsta(iev)
               sumsq = sumsq + d(ista+naddat)**2
             enddo
c  number of degrees of freedom based on amplitude misfit - only two
c  adjustable amplitude parameters per event
c             sigma2 = sumsq/(nsta(iev)-2)
c  number of degrees of freedom assuming 6 wave parameters per event and that
c  each event supplies average of about 2 pieces of information about the 
c  velocity structure
             sigma2 = sumsq/(2*nsta(iev)-8)
             stddevdata(iev) = stddevdata(iev)*sqrt(sigma2)
c  set floor on how small stddevdata is allowed to get - below .03 probably
c  unrealistic even for best data
             if (stddevdata(iev).lt.minstd) stddevdata(iev)=minstd
             naddat = naddat + 2*nsta(iev)
           enddo
           icnt = 1
           go to 100
         endif



c
c  By the time this converges to a solution, d should contain the normalized 
c  misfits
c  Calculate sigma**2 and standard deviation of model parameters
c  Problem with this again with too many parameters - calculate rms instead
c  or use sigma**2 = 1

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ calculate the uncertainty and write out the results $$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


        sumsq = 0.0
        naddat = 0
        totsumsqph  = 0.0
      totsumsqamp = 0.0
        do iev = 1, nevents
          write(10,*) 'event ', idnum(iev)
          sumsqph  = 0.0
        sumsqamp =  0.
        sumsqtemp = 0.0
          do ista = 1,nsta(iev)
            dresid1=d(ista+naddat)*stddevdata(iev)
            dresid2=d(ista+naddat+nsta(iev))*stddevdata(iev)
c  convert misfits back to original amplitude and phase
            predamp = amprms(iev,ifreq)
     1        *sqrt((streal(iev,ista,ifreq)-dresid1)**2 
     2                + (stimag(iev,ista,ifreq)-dresid2)**2)
            prefase = atan2(-(stimag(iev,ista,ifreq)-dresid2),
     1                (streal(iev,ista,ifreq)-dresid1))/twopi
c check for prefase being off by twopi
           if (prefase - dph(iev,ista).gt.0.5) prefase = prefase-1.0
           if (prefase - dph(iev,ista).lt.-0.5) prefase = prefase+1.0   

           write(10,105)  fn(iev,ista),ista,
     2         staamp(iev,ista,ifreq)/amprms(iev,ifreq), 
     2         predamp/amprms(iev,ifreq), predamp/amprms(iev,ifreq) -
     2         staamp(iev,ista,ifreq)/amprms(iev,ifreq),
     2       dph(iev,ista), prefase, prefase-dph(iev,ista),
     2    (prefase-dph(iev,ista))*(1./freq(1))

 105      format(a30,i4, 2x, 6f8.3,f6.1 )

           sumsq = sumsq+(d(ista+naddat)**2+d(ista+naddat+nsta(iev))**2)
     1               *stddevdata(iev)**2
           sumsqtemp = sumsqtemp+(d(ista+naddat)**2+
     1                 d(ista+naddat+nsta(iev))**2)
     1               *stddevdata(iev)**2

           sumsqph = sumsqph + (dph(iev,ista)-prefase)**2
         sumsqamp = sumsqamp + ((predamp-staamp(iev,ista,ifreq))
     1                           /amprms(iev,ifreq))**2

          enddo
          totsumsqph = totsumsqph + sumsqph/freq(ifreq)**2
        totsumsqamp = totsumsqamp + sumsqamp
          rmsphase(iev) = sqrt(sumsqph/nsta(iev))/freq(ifreq)
        rmsamp(iev)   = sqrt(sumsqamp/nsta(iev))
        rmsdata(iev) = sqrt(sumsqtemp/(2*nsta(iev)))
          naddat = naddat + 2*nsta(iev)
        enddo
c  rmsph multiplied by 2 because only half the number of observations
        rmsph = sqrt(2.*totsumsqph/nobs)
      rmsamplitude = sqrt(2.*totsumsqamp/nobs)
        sigma2 = sqrt(sumsq/nobs)
        do j = 1, np
          stddev(j) = sqrt(gtginv(j,j))
        enddo
c        write(10,*) 'nobs',nobs, 'rank', rank, ' 
c     1      rank from vel params ', rank2
c        write(10,*) icnttot, icnttot2,icnt, ' iterations', sigma2, 
c     1     ' unnormalized rms misfit', rmsph, 'rms phase misfit,  s',
c     1        '  rms amp mistfit   ', rmsamplitude

        write(11,*) 'nobs',nobs, 'rank', rank, 
     1      ' rank from vel params ', rank2
        write(11,*) icnttot, icnttot2,icnt, ' iterations', sigma2, 
     1     ' unnormalized rms misfit', rmsph, 'rms phase misfit,  s',
     1        '  rms amp mistfit   ', rmsamplitude

c find median event rms phase misfit
        imed1 = (nevents+1)/2
        imed2 = (nevents+2)/2
        do iev = 1, nevents
           sortrms(iev) = rmsphase(iev)
        enddo
        call shell(nevents, sortrms)
        rmsmedian = (sortrms(imed1)+sortrms(imed2))/2.0
c        write(10,*) 'median event misfit', rmsmedian
        write(11,*) 'median event misfit', rmsmedian
        do iev = 1, nevents
          ip = (iev-1)*6

c          write(10,*) 'event ', idnum(iev), rmsdata(iev), 
c     1  ' data std dev', rmsphase(iev), 'rms phase misfit  s',
c     1  '  rms amp mistfit   ',  rmsamp(iev)

          write(11,*) 'event ', idnum(iev), rmsdata(iev), 
     1  ' data std dev', rmsphase(iev), 'rms phase misfit  s',
     1  '  rms amp mistfit   ',  rmsamp(iev)

          wvaz1= stazim1(iev)/convdeg
          wvaz2 = stazim2(iev)/convdeg
          stdwvaz1 = stddev(ip+3)/convdeg
          stdwvaz2 = stddev(ip+4)/convdeg
c          write(10,*) wvaz1, startamp1(iev), stphase1(iev)
c          write(10,*) wvaz2,startamp2(iev), stphase2(iev)
          write(11,*) wvaz1, startamp1(iev), stphase1(iev)
          write(11,*) wvaz2,startamp2(iev), stphase2(iev)
        enddo

c   write parameter covariance matrix for isotropic velocity parameters
c   for use in calculating variance of velocity model.  Could be smart and
c   write only upper or lower triangle of symmetric matrix to save space.
        open(26,file = fvelarea)
        write(16,*) nnodes
      do ii = 1, nnodes
          do jj = 1, nnodes
          write(16,'(f10.4)') gtginv(ii+i6,jj+i6)
          enddo
        enddo

        do ii = 1, nnodes
c          write(10,*) ii, crrntmod(ii+i6), stddev(ii+i6)
          write(11,*) ii, crrntmod(ii+i6), stddev(ii+i6)
          write(26,*) ii, crrntmod(ii+i6), stddev(ii+i6)
        enddo


      do ii = 1, ntype
       write(11,167)  sttypename(ii), ampmult(ii), phcor(ii) 
      enddo
 167    format(i6, 2f10.3)






      enddo
900   close(unit = 10)
      close(unit = 11)
      close(unit = 12)
      close(unit = 13)
      close(unit = 14)
      close(unit = 15)
      close(unit = 16) 
      close(unit = 26)
      stop 
      end


      SUBROUTINE dlubksb(gtg,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION gtg(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
c      common /gmatrix/ gtg,gtginv

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-gtg(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-gtg(i,j)*b(j)
13      continue
        b(i)=sum/gtg(i,i)
14    continue
      return
      END


      SUBROUTINE dludcmp(gtg,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,gtg(np,np),TINY
      PARAMETER (NMAX=7000,TINY=1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
       
c      common /gmatrix/ gtg,gtginv
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(gtg(i,j)).gt.aamax) aamax=abs(gtg(i,j))
11      continue
c        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=gtg(i,j)
          do 13 k=1,i-1
            sum=sum-gtg(i,k)*gtg(k,j)
13        continue
          gtg(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=gtg(i,j)
          do 15 k=1,j-1
            sum=sum-gtg(i,k)*gtg(k,j)
15        continue
          gtg(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=gtg(imax,k)
            gtg(imax,k)=gtg(j,k)
            gtg(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(gtg(j,j).eq.0.)gtg(j,j)=TINY
        if(j.ne.n)then
          dum=1./gtg(j,j)
          do 18 i=j+1,n
            gtg(i,j)=gtg(i,j)*dum
18        continue
        endif
19    continue
      return
      END

      subroutine disthead(slat,slon,flat,flon,delta,azim)
c  Calculates distance and azimuth on sphere from starting point s 
c  to finishing point f
      dtor= 3.1415928/180.
      slt = slat*dtor
      sln = slon*dtor
      flt = flat*dtor
      fln = flon*dtor
      delta = acos(sin(slt)*sin(flt)+cos(slt)*cos(flt)*cos(fln-sln))
      azim = atan2(sin(fln-sln)*cos(flt),
     1  sin(flt)*cos(slt) - cos(fln-sln)*cos(flt)*sin(slt))
      delta = delta/dtor
      azim = azim/dtor
      return
      end

      subroutine gohead(slat,slon,delta,azim,flat,flon)
c Calculates final latitude and longitude f when starting at point s
c traveling delta degrees on sphere at initial heading azim
      dtor= 3.1415928/180.
      slt = slat*dtor
      dlta = delta*dtor
      azm = azim*dtor
      flat = asin(cos(slt)*sin(dlta)*cos(azm) + sin(slt)*cos(dlta))
      flon = atan2(sin(dlta)*sin(azm), 
     1   cos(slt)*cos(dlta) - sin(slt)*sin(dlta)*cos(azm))
      flat = flat/dtor
      flon = slon + flon/dtor
      if (flon.gt.360.) flon = flon - 360.
      if (flon.lt.-360.) flon = flon + 360.
      return
      end   

      subroutine shell(n,a)
c  shell sort from Press et al.  Sorts a, replacing original, in ascending order
      integer n, i,j, inc
      real a(n), v
      inc = 1
1     inc = 3*inc + 1
      if (inc.le.n) go to 1
2       continue
        inc = inc/3
        do i = inc+1, n
          v = a(i)
          j = i
3         if (a(j-inc).gt.v) then
            a(j) = a(j-inc)
            j = j-inc
            if (j.le.inc) go to 4
          go to 3
          endif
4         a(j) = v
        enddo
        if(inc.gt.1) go to 2
        return
        end




      real function misfit(p)
c   calculates misfit for single event - use for 2-plane wave minimazation
c nnodes

       parameter (maxnfreq=1, maxnsta=120, maxnstatype = 1500,
     1    maxpts=20000, nparam = 5000,
     1    maxnobs =50000, maxnodes = 2000, maxevnts = 500,
     1    maxnxints = 501, maxnyints = 501,ndeg = 81 )

      real*4 streal(maxevnts,maxnsta,maxnfreq)
      real*4 stimag(maxevnts,maxnsta,maxnfreq)
      real*4 rloc(maxevnts,maxnsta), azloc(maxevnts,maxnsta)
      real*4 d(maxnobs), xsta(maxevnts,maxnsta),ysta(maxevnts,maxnsta)
      real*4 freq(maxnfreq)
      real*4 stddevdata(maxevnts)
      real*4 dtime(maxnsta), avslow(maxnsta)
      integer*4 nsta(maxevnts), iref(maxevnts)

      real*4 xbox(maxevnts,4),ybox(maxevnts,4)

      real*4 dxnode,dynode
      real*4 sensitivity(maxnxints,maxnyints)
      real*4 ampsens(maxnxints,maxnyints)
      real*4  damp1per(maxnsta), damp2per(maxnsta)
      real*4  dphase1(maxnsta),dphase2(maxnsta)
      real*4   phase1(maxnsta), phase2(maxnsta)

      real*4 phase(maxnsta,ndeg),dphase(maxnsta,ndeg),
     1      dampper(maxnsta,ndeg)

      real*4  appvel(maxnodes)

      real*4   ampmult(maxnstatype),phcor(maxnstatype)
      integer*4 iptype(maxnstatype)
     
      real*4 xmin(maxevnts,ndeg)


      integer*4 istanum(maxevnts,maxnsta),istacor(maxnsta)
      integer*4   istavar(maxnsta)
      integer*4 istatype(maxnsta)
      integer*4 nxkern
      real*4   xbegkern,dxkern
      integer*4 nykern
      real*4   ybegkern,dykern

      real*4 xnode(maxevnts,maxnodes),ynode(maxevnts,maxnodes)

      real*4 dtime1(maxnsta),dtime2(maxnsta)

      dimension p(6)

      common /residua/ d,rloc,azloc,freq,xsta,dtime,avslow,
     1     streal,stimag,stddevdata,iref,nsta,iev,naddat,ifreq,
     1    xbox,ybox,ysta,nxkern,nykern,dxkern,dykern,dxnode,dynode,
     1   sensitivity,ampsens,unifvel,appvel,ampmult,istavar,istanum,
     1  xnode, ynode,nnodes,xmin,
     1   phase,dphase,dampper,phcor,ntype,iptype,istatype

        

        twopi = 3.1415928*2.
        onepi = 3.1415928
        convdeg = 3.1415928/180.
      misfit = 0.0
      startamp1 = p(1)
      startamp2 = p(2)
        stazim1 = p(3)
        stazim2 = p(4)
        stphase1 = p(5)
        stphase2 = p(6)



       if(stazim1 .lt. -38*convdeg) stazim1 = -38*convdeg
       if(stazim1 .gt. 38*convdeg)  stazim1 = 38*convdeg
       if(stazim2 .lt. -38*convdeg) stazim2 = -38*convdeg
       if(stazim2 .gt. 38*convdeg)  stazim2 = 38*convdeg

         ideg1 = int(stazim1/convdeg+(ndeg-1)/2) + 1
         ideg2 = int(stazim2/convdeg+(ndeg-1)/2) + 1


c  find minimum xbox

       do 10 ista = 1, nsta(iev)


       xstatemp = 
     1   xsta(iev,ista)*cos(stazim1)
     1 + ysta(iev,ista)*sin(stazim1)               
     
       ystatemp = 
     1  -xsta(iev,ista)*sin(stazim1)
     1 + ysta(iev,ista)*cos(stazim1)     

   
c  normalize weights to be equivalent to distance
        
      dtime1(ista) = (xstatemp-xmin(iev,ideg1))*(1.0/unifvel)

  10    enddo


c       calculate for second plane wave

  
        do 20 ista = 1, nsta(iev)
       xstatemp = 
     1   xsta(iev,ista)*cos(stazim2)
     1 + ysta(iev,ista)*sin(stazim2)               
     
       ystatemp = 
     1  -xsta(iev,ista)*sin(stazim2)
     1 + ysta(iev,ista)*cos(stazim2)     

      dtime2(ista) = (xstatemp-xmin(iev,ideg2))*(1.0/unifvel)

  20    enddo

      do ista = 1, nsta(iev)

         phase1(ista) = dtime1(ista)*freq(1)
         phase2(ista) = dtime2(ista)*freq(1)

      enddo



      do ista = 1,nsta(iev)

         iptptemp = iptype(istatype(istanum(iev,ista)))
         ampadj =  ampmult(iptptemp)
          phadj =    phcor(iptptemp)

        prefase1 = ( (phase1(ista)+dphase(ista,ideg1))
     1   - (phase1(iref(iev))+dphase(iref(iev),ideg1)))
     1             + stphase1
        prefase2 = ( (phase2(ista)+dphase(ista,ideg2))
     1   - (phase2(iref(iev))+dphase(iref(iev),ideg2)))
     1             + stphase2

          prefase1=prefase1+phadj
          prefase2=prefase2+phadj


        staamp1 = startamp1*(1.+dampper(ista,ideg1))
        staamp2 = startamp2*(1.+dampper(ista,ideg2))


          cosph1 = cos(prefase1*twopi)
          cosph2 = cos(prefase2*twopi)
          sinph1 = sin(prefase1*twopi)
          sinph2 = sin(prefase2*twopi)
          
          prereal = staamp1*cosph1+staamp2*cosph2  
        preimag = -1.0*(staamp1*sinph1+ staamp2*sinph2)

          prereal = prereal*ampadj
          preimag = preimag*ampadj

c
c  data vector and partial derivatives listed event by event with all
c  real data for first event followed by imaginary data, then onto next event
c  d contains misfit to starting model

c
          kreal = ista + naddat
          kimag = kreal + nsta(iev)
          d(kreal) = (streal(iev,ista,ifreq) - prereal)
     1       /stddevdata(iev)
          d(kimag) = (stimag(iev,ista,ifreq) - preimag)
     2        /stddevdata(iev)
          misfit = misfit + d(kreal)**2 + d(kimag)**2

       enddo
      return 
        end

      integer function blank(file)
      character file*70
      do 50 i=1,70
      if(file(i:i).ne.' ') goto 50
      blank=i-1
      return
50     continue
      write(1,100) file
100   format(' no blanks found in ',a70)
      blank = 0
      return
      end


      integer function findid(inpstaname)
      parameter ( nsta = 112)
      character *6 inpstaname,staname 

      open(22,file = 'stationid.dat',status = 'old')

        findid = 0

 19     read(22,*,end = 111)  staname, idsta
      if ( inpstaname .eq. staname) then
            findid = idsta
            goto 111
      endif             
        goto 19
111     continue 

      if( findid ==0 ) then
         write(*,*) staname, inpstaname,"  ",'station not in the list'
         stop 
      endif


      close(22)
      return
      end   



      integer function findstatype(stanumber)
      parameter ( fnsta = 100)
      integer statype,stanumber
      character*6 dummy
      open(23,file = 'stationid.dat',status = 'old')

      findstatype = 0 
 155    read(23,*,end=114) dummy,idsta,statype
      if ( stanumber .eq. idsta ) then
            findstatype = statype
              goto 114
      endif
      goto 155
 114    continue


      close(23)
      return
      end   


