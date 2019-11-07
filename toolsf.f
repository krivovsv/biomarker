! subroutines to determine optimal biomarker from longitudinal data
! la,lrc,dx2,info=evec(lx,lp) - approximates the second eigenvector by a linear combination of variables - unsupervised determination of optimal biomarker 
! la,lrc,dx2,info=committor(lx,lp,bval,skip) -. determines the commitor function by a linear combination of variables - supervised determination of optimal biomarker
! lrc,dx2,info=npevec(lx,lp) - finds an approximation to the second eigenvector by non-parametric optimization of variables - nonparametric unsupervised determination of optimal biomarker
! lrc,dx2,info=npcommittor(lx,lp,bval) - finds an approximation to the committor function by non-parametric optimization of variables - nonparametric upervised determination of optimal biomarker





! INPUT
! np -number of patients, nt - the total number of samples, summed over all patients, nx -number of variables in a sample
! lx(nx,nt) - variables in the samples, input
! lp(nt) - patient to whom the sample belongs, should look like 1,1,1, ...,1,2,2..
! bval(nt) - bounary value for the sample, to perform the supervised biomarker discovery. it is 0 for all samples but the last sample of a patient. 1 it means that the pation belongs to first group, -1 patient belongs to the second group, 2 - no information about the group.
! skip(nt) - if nonzero, then point is skipped from optimization, while the rc is still computed for the point. used for cross-validations.


! OUTPUT
! la(nx+1) - coefficients of the linear combination of the variables, the last coefficient, la(nx+1) is that before 1.
! lrc(nt) - values of the optimal biomarker for each sample

      subroutine committor(lx,lp,bval,skip,la,lrc,dx2,lmb,info,nx,nt)
! la,lrc,dx2,info=committor(lx,lp,bval) -. determines the commitor function by a linear combination of variables - supervised determination of optimal biomarker
      implicit none
      integer nt,nx,info
      real*8 lx(nx,nt),lrc(nt),dx2,la(nx+1),lmb
      integer lp(nt),bval(nt),skip(nt)

      integer i,iset
      real*8 compdx2
      real*8 y(nx+1,nt)
      
      
Cf2py intent(out):: la,lrc,dx2,info  

      do iset=1,nt
        do i=1,nx
          y(i,iset)=lx(i,iset)
        enddo
        y(nx+1,iset)=1
      enddo        
      call optimdx2lc(y,bval,skip,nx+1,nt,la,lrc,lmb,info)
      dx2=compdx2(lrc,bval,nt)
      end

      subroutine evecs(lx,lp,lrc,skip,info,nx,nt,lmb)
! lrc,info=evecs(lx,lp) - approximates eigenvectors higher than second by a linear combination of variables - unsupervised determination of optimal biomarker 
      implicit none
      integer nt,nx,info
      real*8 lx(nx,nt),lrc(nx,nt),dx2,lmb
      integer lp(nt),skip(nt)

      integer i,iset
      real*8 compdx2evec
      real*8 y(nx+1,nt)
      
      
Cf2py intent(out):: lrc,dx2,info  
      
      do iset=1,nt
        do i=1,nx
          y(i,iset)=lx(i,iset)
        enddo
        y(nx+1,iset)=1
      enddo
      call optimeveclc(y,lp,skip,nx+1,nt,lrc,lmb,info)
      end

      subroutine npcommittor(lx,lp,bval,lrc,dx2,info,nx,nt,np,eps,miter)
! lrc,dx2,info=npcommittor(lx,lp,bval,...) - finds an approximation to the committor function by non-parametric optimization of variables - nonparametric upervised determination of optimal biomarker
! np - the degree of the polynomial
! eps - accuracy of optimization
! miter - maximal number of iterations
      implicit none
      integer nt,nx,info,np,miter
      real*8 lx(nx,nt),lrc(nt),dx2,eps
      integer lp(nt),bval(nt)

      integer iter,i,j
      real*8 y(nt)
      real*8 compdx2,dx2last
      
Cf2py intent(out):: lrc,dx2,info

      lrc=0
      dx2=compdx2(lrc,bval,nt)
      dx2last=dx2
      do iter=1,miter
        i=rand()*nx+1 
        i=min(i,nx)
        do j=1,nt
          y(j)=lx(i,j)
        enddo
        if (iter==1)then
          call optimdx2np(1,1,lrc,y,bval,nt,info,dx2)
        else
          call optimdx2np(np,np,lrc,y,bval,nt,info,dx2)
        endif
c        if (info/=0)return
        if (mod(iter,10)==0 .and. info==0)then
          if (abs(dx2last-dx2)<eps) return
          dx2last=dx2
        endif  
        if (mod(iter,5)==0)call optimdx2np(6,0,lrc,y,bval,nt,info,dx2)
      enddo  
      end

      subroutine optimdx2lc(y,bval,skip,ny,nt,la,rc,lmb,info)
!     optimal linear combination of yi  
!!!!  minimum of [x(t)-x(t+dt)]^2 +[x(t)-x_A]^2+[x(t)-x_B]^2 where x_A=-1 and x_B=1
      implicit none
      integer ny,nt,bval(nt),skip(nt)
      real*8 rc(nt),y(ny,nt),la(ny),lmb

      
      real*8 drdr(ny, ny),dind
      real*8 al(ny),alal,rhs(ny)
      
     
      integer i,j,i1,iset,k,isize
      real*8 r
      integer info,ipiv(ny)
      
      
      integer i12(ny),i21(ny),nij2
      real*8 dfij(ny)
      
      isize=ny
      do i=1,isize
        al(i)=0
        do j=1,i
          drdr(i,j)=0
        enddo
      enddo
      
      alal=0
      do iset=1,nt 
        if (skip(iset)/=0)cycle ! do not consider points
        if (bval(iset)==0)then  ! not a boundary point
          do i=1,ny
            dfij(i)=y(i,iset+1)-y(i,iset)
          enddo  
          do j=1,ny
            do k=1,j
              drdr(j,k)=drdr(j,k)+dfij(k)*dfij(j)
            enddo
          enddo
        else
          if (bval(iset)/=2)then ! a boundary point connected to +1 or -1
            dind=-bval(iset) 
            do i=1,ny
              dfij(i)=y(i,iset)
            enddo  
            alal=alal+dind*dind
            do j=1,ny
              al(j)=al(j)+2*dind*dfij(j)
            enddo
            do j=1,ny
              do k=1,j
                drdr(j,k)=drdr(j,k)+dfij(k)*dfij(j)
              enddo
            enddo
          endif
        endif  
      enddo
      nij2=0
      do i=1,isize
        i21(i)=0
      enddo  
      do i=1,isize  ! select only non zero
        if (drdr(i,i)<1e-7)cycle
        nij2=nij2+1
        i12(nij2)=i
        i21(i)=nij2
        drdr(i,i)=drdr(i,i)+lmb
      enddo
      do i=1,nij2
        do j=1,i
          drdr(i,j)=drdr(i12(i),i12(j))
          drdr(j,i)=drdr(i,j)
        enddo
        al(i)=al(i12(i))
        rhs(i)=-al(i)/2
      enddo  
      call dgetrf(nij2,nij2,drdr,isize,ipiv,info)
      if(info/=0)then
        write(*,*)'dgetrf info=',info
        return
        stop'info/=0'
      endif
      call dgetrs('N',nij2,1,drdr,isize,ipiv,rhs,isize,info)
      if(info/=0)then
        write(*,*)'dgetrs info=',info
        return
        stop'info/=0'
      endif
      
      do iset=1,nt
          do i=1,ny
            dfij(i)=y(i,iset)
          enddo
          r=0
          do i=1,nij2
            i1=i12(i)
            r=r+rhs(i)*dfij(i1)
          enddo
c          if (r>1)r=1     ! impose boundaries 0<q<1
c          if (r<-1)r=-1    
          rc(iset)=r
      enddo
      la=0
      do i=1,nij2
        la(i12(i))=rhs(i) ! assign coefficients
      enddo
      end    
      
      subroutine optimeveclc(y,lp,skip,ny,nt,rc,lmb,info)
!     optimal linear combination of yi      
!!!!  minimum of [x(t)-x(t+dt)]^2 under constraint sum_t x^2(t)=1
      implicit none
      integer ny,nt,lp(nt),skip(nt)
      real*8 rc(ny-1,nt),y(ny,nt),lmb

      
      real*8 drdr(ny, ny),rr(ny,ny)
     
      integer i,j,iset,k
      real*8 r
      integer info,lwork
      real*8 work(ny*10),w(ny)
      
      integer nij2,i12(ny),i21(ny),iev
      real*8 dfij(ny)
      
      drdr=0
      rr=0
      lwork=ny*10
      
      do iset=1,nt
          if (skip(iset)/=0)cycle ! do not consider points
          do i=1,ny
            dfij(i)=y(i,iset)
          enddo  
          do j=1,ny
            do k=1,j
              rr(j,k)=rr(j,k)+dfij(k)*dfij(j)
            enddo
          enddo
      enddo
      do iset=1,nt-1
        if (skip(iset)/=0)cycle ! do not consider points
        if (lp(iset)==lp(iset+1))then  ! the same patient
          do i=1,ny
            dfij(i)=y(i,iset+1)-y(i,iset)
          enddo  
          do j=1,ny
            do k=1,j
              drdr(j,k)=drdr(j,k)+dfij(k)*dfij(j)
            enddo
          enddo
        endif
      enddo
      nij2=0
      i21=0
      do i=1,ny  ! select only non zero
        if (drdr(i,i)<1e-5 .and. i<ny)cycle
        nij2=nij2+1
        i12(nij2)=i
        i21(i)=nij2
        drdr(i,i)=drdr(i,i)+lmb
        rr(i,i)=rr(i,i)+lmb
      enddo
      do i=1,nij2
        do j=1,i
          drdr(i,j)=drdr(i12(i),i12(j))
          drdr(j,i)=drdr(i,j)
          rr(i,j)=rr(i12(i),i12(j))
          rr(j,i)=rr(i,j)
        enddo
      enddo
      call dsygv(1,'V','U',nij2,drdr,ny,rr,ny,w,work,lwork,info)
      if (info/=0)then
        write(*,*)'info=',info
        write(*,*)'ny,nij2=',ny,nij2
        return
      endif
      do i=1,ny
        write(*,*)i,w(i)
      enddo
      do iev=1,ny-1
        do iset=1,nt
          r=0
          do i=1,nij2
            r=r+y(i12(i),iset)*drdr(i,iev)
          enddo
          rc(iev,iset)=r
        enddo
      enddo
      end    

      
      real*8 function compdx2(rc,rcbd,nsets)
      implicit none
      integer nsets,rcbd(nsets)
      real*8 rc(nsets)
      
      integer iset,dt
      real*8 r
      
      r=0
      dt=1
      do iset=1,nsets
           if (rcbd(iset)==0)then
            r=r+(rc(iset+dt)-rc(iset))**2
          else
            if (rcbd(iset)/=2) r=r+(rcbd(iset)-rc(iset))**2
          endif  
      enddo
      compdx2=r
      end
      
      real*8 function compdx2evec(rc,lp,nsets)
      implicit none
      integer nsets,lp(nsets)
      real*8 rc(nsets)
      
      integer iset
      real*8 r
      
      r=0
      do iset=1,nsets-1
           if (lp(iset)==lp(iset+1))then
            r=r+(rc(iset+1)-rc(iset))**2
          endif  
      enddo
      compdx2evec=r
      end


      subroutine optimdx2np(nr,ny,rc,y,bval,nt,info,dx2)
      implicit none
      integer nr,ny,nt,info,bval(nt)
      real*8 rc(nt),y(nt),dx2

!!!!  minimum of [x(t)-x(t+dt)]^2 +[x(t)-x_A]^2+[x(t)-x_B]^2 where x_A=-1 and x_B=1
!!!!  for all i: \sum_j al_j[<dx_idx_j>_t+<x_ix_j>_AB] = x_AB <x_i>_AB 
      
      real*8 dxdx((ny+1)*(ny+1)+nr-ny,(ny+1)*(ny+1)+nr-ny)
      real*8 al((ny+1)*(ny+1)+nr-ny),alal
      real*8 rhs((ny+1)*(ny+1)+nr-ny)
      
     
      integer i,j,i1,nij,iset,isize
      real*8 r
      
      integer ipiv((ny+1)*(ny+1)+nr-ny)
      real*8 sval((ny+1)*(ny+1)+nr-ny),work(((ny+1)*(ny+1)+nr-ny)*10)
      integer rank,lwork,iwork(((ny+1)*(ny+1)+nr-ny)*10)
      
      integer i12((ny+1)*(ny+1)+nr-ny)
      integer i21((ny+1)*(ny+1)+nr-ny),nij2
      
      real*8 compdx2,dx2n
      real*8 rcn(nt)

      real*8 rcp1,yp1,dfij((ny+1)*(ny+1)+nr-ny)
     
      isize=(ny+1)*(ny+1)+nr-ny
      lwork=isize*10
      al=0
      alal=0
      dxdx=0
      
      call compxxnp(nr,ny,rc,y,bval,nt,alal,al,dxdx,isize)
      nij2=0
      do i=1,isize
        i21(i)=0
      enddo  
      do i=1,isize  ! select only non zero
        if (dxdx(i,i)<1e-5)cycle
        nij2=nij2+1
        i12(nij2)=i
        i21(i)=nij2
      enddo
      do i=1,nij2
        do j=1,i
          dxdx(i,j)=dxdx(i12(i),i12(j))
          dxdx(j,i)=dxdx(i,j)
        enddo
        al(i)=al(i12(i))
        rhs(i)=-al(i)/2
      enddo    


      if (.true.)then !solving SLE
        call dgetrf(nij2,nij2,dxdx,isize,ipiv,info)
        if(info/=0)then
          write(*,*)'dgetrf info=',info
          return
          stop'info/=0'
        endif
        call dgetrs('N',nij2,1,dxdx,isize,ipiv,rhs,isize,info)
        if(info/=0)then
          write(*,*)'dgetrs info=',info
          stop'info/=0'
        endif
      else           ! solving least square problem
        call DGELSD(nij2,nij2,1,dxdx,isize,rhs,isize,sval,-1d-14, rank,
     $                         work,lwork,iwork,info)   
        if(info/=0)then
          write(*,*)'DGELSD info=',info
          stop'info/=0'
        endif
      endif
      

!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,nij2,ny,nr,y,rc,rcfix,rcn,rhs,i12) 
!$OMP&   PRIVATE(nij,rcp1,yp1,dfij,i1,r)

!$OMP DO
      do iset=1,nt
          rcn(iset)=rc(iset)
          nij=0
          rcp1=1
          do i=0,ny
            yp1=rcp1
            do j=0,ny-i
              nij=nij+1
c              dfij(nij)=rc(iset)**i*y(iset)**j
              dfij(nij)=yp1
              yp1=yp1*y(iset)
            enddo
            rcp1=rcp1*rc(iset)
          enddo  
          do i=ny+1,nr
            nij=nij+1
c            dfij(nij)=rc(iset)**i
            dfij(nij)=rcp1
            rcp1=rcp1*rc(iset)
          enddo    
          r=0
          do i=1,nij2
            i1=i12(i)
            r=r+rhs(i)*dfij(i1)
          enddo    
          rcn(iset)=r
      enddo
!$OMP END PARALLEL
      dx2n=compdx2(rcn,bval,nt)
      if (.not. (dx2n>1 .or. dx2n<1)) return ! check for nan
      if (dx2n>dx2) return
      rc=rcn
      dx2=dx2n 
      end    


      subroutine compxxnp(nr,ny,rc,y,bval,nt,alal,al,dxdx,isize2)
      implicit none
      integer nt,isize2,nr,ny,bval(nt)
      real*8 rc(nt),y(nt),alal,al(isize2),dxdx(isize2,isize2)
      
      integer iset,nij,i,j
      real*8 rcp1,rcp2,yp1,yp2,dfij(isize2)
      real*8 dind
      
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,rcind,rcfix,ny,nr,y,rc) 
!$OMP&   PRIVATE(bi,b1,b2,dind,ix,nij,rcp2,rcp1,yp2,yp1,dfij)
!$OMP&   reduction(+:alal,al,dxdx) 

!$OMP DO
        do iset=1,nt
  !       f(t)=2rcind(t)+(-1)**rcind(t)*rc(t)=k+b*rc
          if (bval(iset)==0)then
            nij=0
            rcp2=1
            rcp1=1
            do i=0,ny
              yp2=rcp2
              yp1=rcp1
              do j=0,ny-i
                nij=nij+1
c                dfij(nij)=rc(iset+dt)**i*y(iset+dt)**j
c     $                     -rc(iset)**i*y(iset)**j
                 dfij(nij)=yp2-yp1
                 yp2=yp2*y(iset+1)
                 yp1=yp1*y(iset)
              enddo
              rcp2=rcp2*rc(iset+1)
              rcp1=rcp1*rc(iset)
            enddo
            do i=ny+1,nr
              nij=nij+1
c              dfij(nij)=rc(iset+dt)**i-rc(iset)**i
              dfij(nij)=rcp2-rcp1
              rcp2=rcp2*rc(iset+1)
              rcp1=rcp1*rc(iset)
            enddo   
            do i=1,nij
              do j=1,i
                dxdx(i,j)=dxdx(i,j)+dfij(i)*dfij(j)
              enddo
            enddo
          else
            if (bval(iset)==2)cycle ! no boundary value 
            dind=-bval(iset)    ! distance offset
            nij=0
            rcp1=1
            do i=0,ny
              yp1=rcp1
              do j=0,ny-i
                nij=nij+1
c                dfij(nij)=bi*rc(ix)**i*y(ix)**j
                dfij(nij)=yp1
                yp1=yp1*y(iset)
              enddo
              rcp1=rcp1*rc(iset)
            enddo  
            do i=ny+1,nr
              nij=nij+1
c              dfij(nij)=bi*rc(ix)**i
              dfij(nij)=rcp1
              rcp1=rcp1*rc(iset)
            enddo    
            alal=alal+dind*dind
            do j=1,nij
              al(j)=al(j)+2*dind*dfij(j)
            enddo
            do i=1,nij
              do j=1,i
                dxdx(i,j)=dxdx(i,j)+dfij(i)*dfij(j)
              enddo
            enddo
          endif  
        enddo
!$OMP END PARALLEL
      end

