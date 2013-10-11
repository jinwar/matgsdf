c The subroutine used the simplex method to do simulated annealing
c from Numerical Recipe pp. 445-447
c subroutine called: misfit ran1 and amotsa

      SUBROUTINE anneal(p,y,mp,np,ndim,pb,yb,ftol,iter,temptr)
      integer iter,mp,ndim,np,nmx
      real ftol,temptr,yb,p(mp,np),y(6),pb(1),misfit
      parameter (nmx=200)
      integer i,idum,ihi,ilo,inhi,j,m,n
      real rtol,sum,swap,tt,yhi,ylo,ynhi,ysave,yt,ytry,
     # psum(nmx),amotsa,ran1
      common /ambsa/ tt,idum
      external misfit
       
      tt=-temptr
   1  do n=1,ndim
         sum=0.
         do m=1,ndim+1
            sum=sum+p(m,n)
         end do
         psum(n)=sum
      end do
   2  ilo=1
      inhi=1
      ihi=2
      ylo=y(1)+tt*log(ran1(idum))
      ynhi=ylo
      yhi=y(2)+tt*log(ran1(idum))
      if (ylo.gt.yhi) then
         ihi=1
         inhi=2
         ilo=2
         ynhi=yhi
         yhi=ylo
         ylo=ynhi
      end if

      do i=3,ndim+1
         yt=y(i)+tt*log(ran1(idum))
         if (yt.le.ylo) then
            ilo=i
            ylo=yt
         end if
         if (yt.gt.yhi) then
            inhi=ihi
            ynhi=yhi
            ihi=i
            yhi=yt
         else if (yt.gt.ynhi) then
            inhi=i
            ynhi=yt
         end if
      end do

      rtol=2.*abs(yhi-ylo)/(abs(yhi)+abs(ylo))
c      write (*,'(''rtol='',e15.6)')rtol
c compute the fractional range from highest to lowest and return if satisfactory
c if returning, put the best pt and value in slot 1
  
      if ((rtol.lt.ftol).or.(iter.lt.0)) then
         swap=y(1)
         y(1)=y(ilo)
         y(ilo)=swap
         do n=1,ndim
            swap=p(1,n)
            p(1,n)=p(ilo,n)
            p(ilo,n)=swap
         end do
         return
       end if
       iter=iter-2
c       write(*,'(''call amotsa and fac=-1'',i3)')iter
       ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,misfit,ihi,yhi,-1.0)
c       if (ier.lt.0) pause 'error occurs at fucntion amotsa'
       if (ytry.le.ylo) then
c       write(*,*)'call amotsa and fac=2.'
          ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,misfit,ihi,yhi,2.0)
c          if (ier.lt.0) go to 1
       else if (ytry.ge.ynhi) then
          ysave=yhi
c          write(*,'(''call amotsa and fac=.5'',i3)')iter
          ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,misfit,ihi,yhi,0.5)
c          if (ier.lt.0) go to 1
          if (ytry.ge.ysave) then
             do i=1,ndim+1
                if (i.ne.ilo) then
                   do j=1,ndim
                      psum(j)=0.5*(p(i,j)+p(ilo,j))
                      p(i,j)=psum(j)
                   end do
                   y(i)=misfit(psum,0)
                end if
             end do
           iter=iter-ndim
           go to 1
          end if
       else
          iter=iter+1
       end if
       go to 2
       END

      FUNCTION amotsa(p,y,psum,mp,np,ndim,pb,yb,misfit,ihi,yhi,fac)
      common /ambsa/ tt,idum
      integer ihi,mp,ndim,np,nmx
      real amotsa,fac,yb,yhi,p(mp,1),pb(1),psum(1),y(6),misfit
      parameter (nmx=200)
      external misfit
      integer idum,j
      real fac1,fac2,tt,yflu,ytry,ptry(nmx),ran1

       fac1=(1.-fac)/ndim
       fac2=fac1-fac
       do j=1,ndim
          ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
       end do
       ytry=misfit(ptry,0)
       if (ytry.le.yb) then
          do j=1,ndim
             pb(j)=ptry(j)
          end do
          yb=ytry
       end if
       yflu=ytry-tt*log(ran1(idum))
       if (yflu.lt.yhi) then
          y(ihi)=ytry
          yhi=yflu
          do j=1,ndim
             psum(j)=psum(j)-p(ihi,j)+ptry(j)
             p(ihi,j)=ptry(j)
          end do
       end if
       amotsa=yflu
       return
       END 
