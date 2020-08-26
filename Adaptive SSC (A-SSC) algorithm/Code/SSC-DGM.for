c=============================================================
c DGM based incremental semi-supervised clustering algorithm
c=============================================================
c     main program
      PARAMETER(maxvar=7000, maxrec=100000, maxclust=200, maxnft=700
     1 ,maxclass=100, maxcannot=10000, maxmust=10000, maxsize=7000000)
      implicit double precision(a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1,x2(maxsize),amed(maxnft),plabel(maxclass),xbest(maxvar),z(maxvar)
     2 ,fval(maxrec),x5(maxvar),x3(maxvar)
      integer nob(maxclass),lcan(maxcannot,2),lmust(maxmust,2)
     1 ,listcl(maxrec),listml(maxrec),listbar(maxrec)
      common /csize/m,/c22/a,/anclust/nclust,/cnft/nft,/cmf/mf
     1 ,/crecord/nrecord,/cnc/nc,/adminim/dminim,/cns/ns
     2 ,/cnob/plabel,nob,/ctnorm/tnorm,/cgamma/gamma1,gamma2
     3 ,/cnpnout/npurity,/cnclass/nclass,/cncan/ncan,lcan
     4 ,/cnmust/nmust,lmust,/cnbar/nbar,listbar
      open(39,file='Results1.txt')
      open(42,file='Results2.txt')
      open(78,file='datainput.txt',status='old',form='formatted')
      open(79,file='cannotlink.txt',status='old',form='formatted')
      open(80,file='mustlink.txt',status='old',form='formatted')
c========================================================================
      PRINT *,' '
      PRINT *,'Number of features:'
      read *,nft
      PRINT *,' '
      PRINT *,'Class outputs?'
      PRINT *,' '
      PRINT *,'   1 - if yes'
      PRINT *,'   2 - if no'
      read *,npurity
      PRINT *,' '
      if (npurity.eq.1) then
         PRINT *,'Number of classes in your dataset:'
         read *,nclass
         PRINT *,' '
         PRINT *,'Output column:'
         read *,noutcom
         PRINT *,' '
      end if
      PRINT *,'Maximum number of clusters:'
      read *,nclust
      PRINT *,' '
c===========================================================
      tlimit=7.2d+04
      do i=1,maxrec
       read(78,*,END=901,ERR=900) (a(i,k),k=1,nft)
       nrecord=i
      end do
  900 stop 'Error in input file'       
  901 WRITE(39,*) 'Input complete. Number of records: ',nrecord
      WRITE(39,*)
      
      do i=1,maxcannot
       read(79,*,END=903,ERR=902) (lcan(i,k),k=1,2)
       ncan=i
      end do
  902 stop 'Error in input file'       
  903 WRITE(39,*) 'Input complete. Number of cannot-links: ',ncan
      WRITE(39,*)

      do i=1,maxmust
       read(80,*,END=905,ERR=904) (lmust(i,k),k=1,2)
       nmust=i
      end do
  904 stop 'Error in input file'       
  905 WRITE(39,*) 'Input complete. Number of must-links: ',nmust
      WRITE(39,*)

      IF(npurity.eq.1) mf=nft-1
      IF(npurity.eq.2) mf=nft
      IF(npurity.eq.2) noutcom=0
      tnorm=0.0d+00
c==============================================================
      if ((npurity.eq.1).AND.(noutcom.gt.0)) then
       do i=1,nrecord
        j1=0
        do j=1,nft
         if (j.ne.noutcom) then
          j1=j1+1
          amed(j1)=a(i,j)
         end if
        end do
        amed(nft)=a(i,noutcom)
        do j=1,nft
         a(i,j)=amed(j)
        end do
       END do
      end if
c===========================================================
      IF(npurity.eq.1) THEN
       do i=1,nclass
        nob(i)=0
       end do
       a2=a(1,nft)
       plabel(1)=a2
       nob(1)=1
       n2=1
       do i=2,nrecord
        clabel=a(i,nft)
        do j=1,n2
         IF(clabel.eq.plabel(j)) GO TO 22
        end do
        n2=n2+1
        plabel(n2)=clabel
        j=n2
 22     nob(j)=nob(j)+1
       end do
       do i=1,nclass
         WRITE(39,41) i,plabel(i),nob(i)
       end do
      END if
 41   FORMAT('Class',i4,' label',f7.0,'  number of instances',i7)
c======================================================================
       WRITE(39,*)
       WRITE(39,701)
 701   FORMAT('#Clust','            fval   ','         DB',
     1 '     purity','    v_tot','   v_cant','   v_must','      sil',
     2 '     NMI')
       WRITE(39,*)
c==================================================================
c outliers and normals
c===================================================================
      nucl=0
      numl=0
      nbar=0
      do i=1,nrecord
       do j=1,ncan
        if((i.eq.lcan(j,1)).or.(i.eq.lcan(j,2))) then
         nucl=nucl+1
         listcl(nucl)=i
         go to 115
        end if
       end do
       do j=1,nmust
        if((i.eq.lmust(j,1)).or.(i.eq.lmust(j,2))) then
         numl=numl+1
         listml(numl)=i
         go to 115
        end if
       end do
       nbar=nbar+1
       listbar(nbar)=i
 115  end do
c======================================================================
      if(nrecord.le.200) then
       gamma1=1.0d+00
       gamma2=1.0d+00
       gamma3=1.0d+00
      end if

      if((nrecord.gt.200).and.(nrecord.le.1000)) then
       gamma1=8.0d-01
       gamma2=5.0d-01
       gamma3=4.0d-01
      end if
      
      if((nrecord.gt.1000).and.(nrecord.le.5000)) then
       gamma1=7.0d-01
       gamma2=4.0d-01
       gamma3=3.0d-01
      end if

      if((nrecord.gt.5000).and.(nrecord.le.15000)) then
       gamma1=5.0d-01
       gamma2=3.0d-01
       gamma3=2.0d-01
      end if
      if((nrecord.gt.15000).and.(nrecord.le.50000)) then
       gamma1=4.0d-01
       gamma2=2.0d-01
       gamma3=1.0d-01
      end if
      if(nrecord.gt.50000) then
       gamma1=3.0d-01
       gamma2=1.0d-01
       gamma3=1.0d-01
      end if
      call cpu_time(time1)
c==========================================
      do nc=1,nclust
       PRINT 42,nc
 42    FORMAT('Cluster No.:',i10)
       if(nc.eq.1) then
                   call step1(f,x)
                   toler=1.0d-02*f/dble(nrecord)
                   go to 1
       END if
       call step2(toler,nstart,x2)
       m=mf
       fbarmin=1.0d+26
       fbarmax=0.0d+00
       do j=1,nstart
        do k=1,mf
         z(k)=x2(k+(j-1)*mf)
        end do
        ns=1
        call dgm(z,barf)
        fval(j)=barf
        fbarmin=dmin1(fbarmin,barf)
        fbarmax=dmax1(fbarmax,barf)
        do k=1,mf
         x2(k+(j-1)*mf)=z(k)
        end do
       end do
c==============================================
       fbarmin=fbarmin+gamma3*(fbarmax-fbarmin)
       nstart1=0
       do j=1,nstart
        if (fval(j).le.fbarmin) then
         nstart1=nstart1+1
         do k=1,mf
          x5(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
         end do
        end if
       end do

       nstart=nstart1
       do i=1,nstart
        do k=1,mf
         x2(k+(i-1)*mf)=x5(k+(i-1)*mf)
        end do
       end do
c==============================================
       do k=1,mf
        x5(k)=x2(k)
       end do
       nstart2=1
       do j=2,nstart
        do j1=1,nstart2
         f31=0.0d+00
         do k=1,mf
          f31=f31+(x5(k+(j1-1)*mf)-x2(k+(j-1)*mf))**2
         end do
         IF(f31.LE.toler) GO TO 1200
        end do
        nstart2=nstart2+1
        do k=1,mf
         x5(k+(nstart2-1)*mf)=x2(k+(j-1)*mf)
        end do
1200   end do
       do i=1,nstart2
        do k=1,mf
         x2(k+(i-1)*mf)=x5(k+(i-1)*mf)
        end do
       end do
       nstart=nstart2
c==============================================
       m=mf*nc
       fbest=1.0d+28
       do j=1,nstart
        do i=1,mf
         x(i+(nc-1)*mf)=x2(i+(j-1)*mf)
        END do
        do j1=1,m
         x3(j1)=x(j1)
        end do
        ns=2
        call dgm(x3,fcurrent)
        if (fcurrent.lt.fbest) then
         fbest=fcurrent
         do j1=1,m
          xbest(j1)=x3(j1)
         end do
        end if
       end do
       f=fbest
       do j1=1,m
        x(j1)=xbest(j1)
       end do
c================================================================
  1    continue
       if(nc.gt.1) call clusters(x,f)
       if(npurity.eq.2) purity=0.0d+00
       call cpu_time(time2)
       write(42,603) nc,f,tnorm,time3
 603   format(i6,f30.4,f18.0,f11.3)

       call finresult(x)
       time3=time2-time1
       IF(time3.gt.tlimit) GO TO 2
      end do
   2  continue
      close(39)
      close(42)
      close(78)
      close(79)
      close(80)
      stop
      end

      subroutine clusters(x,f)
      PARAMETER(maxdim=7000, maxrec=100000, maxclust=200, maxnft=700,
     1 maxclass=100, maxcannot=10000, maxmust=10000)
      implicit double precision (a-h,o-z)
      double precision x(maxdim),a(maxrec,maxnft),dminim(maxrec)
     1 ,plabel(maxclass),rad(maxclust),radmax(maxclust),ratio(maxclust)
      integer nob(maxclass),lcand(maxrec),list1(maxrec),lcand1(maxrec)
     1 ,lcan(maxcannot,2),lmust(maxmust,2),listbar(maxrec),nel(maxclust)
     2 ,nk(maxclust,maxrec)
      common /c22/a,/cmf/mf,/anclust/nclust,/crecord/nrecord,/cnc/nc
     1 ,/cnob/plabel,nob,/cnclass/nclass,/cnpnout/npurity
     2 ,/ccand/ncand,lcand,/adminim/dminim,/ctnorm/tnorm,/clist1/list1
     3 ,/cncan/ncan,lcan,/cnmust/nmust,lmust,/cnbar/nbar,listbar

      do j=1,nclust
       nel(j)=0
       rad(j)=0.0d+00
       radmax(j)=0.0d+00
      END do

c===================================================
c Set Abar
c===================================================
      f=0.0d+00
      do k3=1,nbar
       k=listbar(k3)
       f2=1.0d+22
       do j=1,nc
        f1=0.0d+00
        do k1=1,mf
         f1=f1+(a(k,k1)-x(k1+(j-1)*mf))**2
        END do
        tnorm=tnorm+1.0d+00
        if(f2.gt.f1) then
         f2=f1
         jmin=j
        end if
       END do
       dminim(k)=f2
       f=f+f2
       list1(k)=jmin
       rad(jmin)=rad(jmin)+f2
       radmax(jmin)=dmax1(radmax(jmin),f2)
      end do
c===================================================
c Set A_cl - cannot links
c===================================================
      do k=1,ncan
       k1=lcan(k,1)
       k2=lcan(k,2)
       f4=1.0d+22
       do j=1,nc
        do i=1,nc
         if(i.ne.j) then
          f1=0.0d+00
          do k4=1,mf
           f1=f1+(a(k1,k4)-x(k4+(j-1)*mf))**2
          END do
          f2=0.0d+00
          do k4=1,mf
           f2=f2+(a(k2,k4)-x(k4+(i-1)*mf))**2
          END do
          tnorm=tnorm+2.0d+00
          f3=f1+f2
          if(f4.gt.f3) then
           f4=f3
           list1(k1)=j
           list1(k2)=i
           dminim(k1)=f1
           dminim(k2)=f2
          end if         
         end if
        end do 
       end do
       f=f+f4
       j1=list1(k1)
       j2=list1(k2)
       rad(j1)=rad(j1)+dminim(k1)
       rad(j2)=rad(j2)+dminim(k2)
       radmax(j1)=dmax1(radmax(j1),dminim(k1)) 
       radmax(j2)=dmax1(radmax(j2),dminim(k2)) 
      end do 
c===================================================
c Set A_ml - must links
c===================================================
      do k=1,nmust
       k1=lmust(k,1)
       k2=lmust(k,2)
       f4=1.0d+22
       do j=1,nc
         f1=0.0d+00
         do k4=1,mf
          f1=f1+(a(k1,k4)-x(k4+(j-1)*mf))**2
         END do
         f2=0.0d+00
         do k4=1,mf
          f2=f2+(a(k2,k4)-x(k4+(j-1)*mf))**2
         END do
         tnorm=tnorm+2.0d+00
         f3=f1+f2
         if(f4.gt.f3) then
          f4=f3
          list1(k1)=j
          list1(k2)=j
          dminim(k1)=f1
          dminim(k2)=f2
         end if         
       end do
       f=f+f4
       j1=list1(k1)
       rad(j1)=rad(j1)+dminim(k1)+dminim(k2)
       radmax(j1)=dmax1(radmax(j1),dminim(k1)) 
       radmax(j1)=dmax1(radmax(j1),dminim(k2)) 
      end do
      if(nc.eq.1) then
       do i=1,nrecord
        list1(i)=1
       end do
      end if   
c===============================================================
      do i=1,nrecord
       k=list1(i)
       nel(k)=nel(k)+1
       nk(k,nel(k))=i
      end do
c===============================================================
      do k=1,nc
       if(nel(k).gt.0) then
        rad(k)=rad(k)/dble(nel(k))
        else 
          rad(k)=0.0d+00
       end if
      end do

      if(nrecord.lt.500) then
       do k=1,nc
        rad(k)=0.0d+00
       end do
      end if

      if((nc.gt.5).and.(nrecord.gt.500)) then
       ratmin=1.0d+26
       do k=1,nc
        ratio(k)=radmax(k)/rad(k)
        ratmin=dmin1(ratmin,ratio(k))
       end do      
       do k=1,nc
        step1=5.0d-01*ratmin/ratio(k)
        rad(k)=rad(k)+step1*(radmax(k)-rad(k))
       end do
      end if 

      if(nc.lt.nclust) then
       ncand=0
       do k=1,nc
        if(nel(k).gt.2) then
         toler3=5.0d-01*rad(k)
         ncand1=0
         do i=1,nel(k)
          i1=nk(k,i)
          if(ncand1.eq.0) then
           if(dminim(i1).gt.rad(k)) then
            ncand=ncand+1
            ncand1=ncand1+1
            lcand(ncand)=i1
            lcand1(ncand1)=i1           
           end if
          end if
          if(ncand1.gt.0) then
           if(dminim(i1).gt.rad(k)) then
            do j=1,ncand1
             j1=lcand1(j)
             d4=0.0d+00
             do k1=1,mf
              d4=d4+(a(i1,k1)-a(j1,k1))**2
             end do
             tnorm=tnorm+1.0d+00
             if(d4.le.toler3) go to 11
            end do
            ncand=ncand+1
            lcand(ncand)=i1
            ncand1=ncand1+1
            lcand1(ncand1)=i1
           end if
          end if
  11     end do
        end if
       end do 
      end if
c============================================================
      return
      end

c=====================================================
c Step1 calculates the center of the dataset
c=====================================================
      subroutine step1(f,x)
      PARAMETER(maxvar=7000, maxrec=100000, maxclust=200, maxnft=700)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft)
      common /c22/a,/cmf/mf,/crecord/nrecord
      do i=1,mf
       x(i)=0.0d+00
       do j=1,nrecord
        x(i)=x(i)+a(j,i)
       END do
       x(i)=x(i)/dble(nrecord)
      END do
      f=0.0d+00
      do i=1,nrecord
       f1=0.0d+00
       do j=1,mf
        f1=f1+(a(i,j)-x(j))*(a(i,j)-x(j))
       END do
       f=f+f1
      END do
      return
      end

c=========================================================================
c  Step2 computes clusters for each data point
c=========================================================================
      subroutine step2(toler,nstart,x2)
      PARAMETER(maxdim=7000, maxrec=100000, maxclust=200, maxnft=700
     1 ,maxsize=7000000)
      implicit double precision (a-h,o-z)
      double precision x2(maxsize),a(maxrec,maxnft),dminim(maxrec)
     1 ,fmin1(maxrec),x4(maxnft),fval(maxrec)
      integer l4(maxrec),lcand(maxrec),lcand1(maxrec)
      common /c22/a,/cmf/mf,/adminim/dminim,/crecord/nrecord
     1 ,/ctnorm/tnorm,/ccand/ncand,lcand,/cgamma/gamma1,gamma2

      if(nrecord.le.200) then
       do i1=1,ncand
        i=lcand(i1)
        do k=1,mf
         x2(k+(i1-1)*mf)=a(i,k)
        end do
       end do
       nstart=ncand
       return
      end if

      nstart=0
      fmin=0.0d+00
      fmax=-1.0d+26
      do i1=1,ncand
       i=lcand(i1)
       d21=0.0d+00
       do l=1,nrecord
        d3=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d3=d3+(a(i,k)-a(l,k))**2
        end do
        d21=d21+dmin1(0.0d+00,d3-dminim(l))
       end do
       fmin1(i)=d21
       fmin=dmin1(fmin,d21)
       fmax=dmax1(fmax,d21)
      end do

      fmin2=fmin+gamma1*(fmax-fmin)
      ncand1=0
      do i1=1,ncand
       i=lcand(i1)
       if (fmin1(i).le.fmin2) then
        ncand1=ncand1+1
        lcand1(ncand1)=i
       end if
      end do

      ncand=ncand1
      do i1=1,ncand
       lcand(i1)=lcand1(i1)
      end do

      do i1=1,ncand
       i=lcand(i1)
       nclose=0
       do j=1,nrecord
        d1=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d1=d1+(a(i,k)-a(j,k))**2
        end do
        IF(d1.LT.dminim(j)) then
         nclose=nclose+1
         l4(nclose)=j
        END if
       end do

       IF(nclose.eq.0) GO TO 1
       do k=1,mf
        d3=0.0d+00
        do j=1,nclose
         j1=l4(j)
         d3=d3+a(j1,k)
        end do
        x4(k)=d3/DBLE(nclose)
       end do
       do j=1,nstart
        d4=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d4=d4+(x2(k+(j-1)*mf)-x4(k))**2
        end do
        IF(d4.le.toler) GO TO 1
       end do
       nstart=nstart+1
       do k=1,mf
        x2(k+(nstart-1)*mf)=x4(k)
       end do
   1  end do
      d2=0.0d+00
      d6=-1.0d+26
      do j=1,nstart
       d21=0.0d+00
       do l=1,nrecord
        d3=0.0d+00
        tnorm=tnorm+1.0d+00
        do k=1,mf
         d3=d3+(x2(k+(j-1)*mf)-a(l,k))**2
        end do
        d21=d21+dmin1(0.0d+00,d3-dminim(l))
       end do
       fval(j)=d21
       d2=dmin1(d2,d21)
       d6=dmax1(d6,d21)
      end do
      
      d2=d2+gamma2*(d6-d2)      

      nstart1=0
      do j=1,nstart
       if (fval(j).le.d2) then
        nstart1=nstart1+1
        do k=1,mf
         x2(k+(nstart1-1)*mf)=x2(k+(j-1)*mf)
        end do
       end if
      end do
      nstart=nstart1
      return
      end

      subroutine auxfunc(f,x)
      PARAMETER(maxvar=7000, maxrec=100000, maxclust=200, maxnft=700
     1 ,maxcannot=10000, maxmust=10000)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,dist1(maxrec) 
      integer iact(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim
     1 ,/ctnorm/tnorm,/cic/iact

       do i=1,nrecord
        dist1(i)=0.0d+00
        do j=1,mf
         dist1(i)=dist1(i)+(x(j)-a(i,j))**2
        end do
       end do
       tnorm=tnorm+dble(nrecord)
       f=0.0d+00
       do i=1,nrecord
        f5=dist1(i)
        f4=dmin1(dminim(i),f5)
        if(f4.eq.dminim(i)) iact(i)=1
        if(f4.eq.f5) iact(i)=2
        f=f+f4
       end do
      return
      end

      subroutine auxgrad(x,grad)
      PARAMETER(maxvar=7000, maxrec=100000, maxclust=200, maxnft=700
     1 ,maxcannot=10000, maxmust=10000)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dminim(maxrec)
     1 ,grad(maxvar)
      integer iact(maxrec)  
      common /c22/a,/cmf/mf,/crecord/nrecord,/adminim/dminim
     1 ,/cic/iact,/cns/ns
     
      do i=1,mf
       grad(i)=0.0d+00
      end do
      do i=1,nrecord
       if(iact(i).eq.2) then
        do j=1,mf
         grad(j)=grad(j)+2.0d+00*(x(j)-a(i,j))
        end do
       end if 
      end do
      return
      end

c=====================================================================
      subroutine func(f,x)
      PARAMETER(maxvar=7000, maxrec=100000, maxclust=200, maxnft=700
     1 ,maxcannot=10000, maxmust=10000)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),dist2(maxrec,maxclust)
      integer lact3(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/ctnorm/tnorm
     2 ,/cic1/lact3
     
      do i=1,nrecord
       do j=1,nc
        dist2(i,j)=0.0d+00
        do k=1,mf
         dist2(i,j)=dist2(i,j)+(x(k+(j-1)*mf)-a(i,k))**2
        end do
       end do
      end do
      tnorm=tnorm+dble(nrecord)*dble(nc)
   
      f=0.0d+00
      do i=1,nrecord
       f4=1.0d+26
       do k=1,nc
        f5=dist2(i,k)
        if(f4.gt.f5) then
         f4=f5
         lact3(i)=k
        end if                 
       end do
       f=f+f4
      end do
      return
      end

      subroutine funcgrad(x,grad)
      PARAMETER(maxvar=7000, maxrec=100000, maxclust=200, maxnft=700
     1 ,maxcannot=10000, maxmust=10000)
      implicit double precision (a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),grad(maxvar)
      integer lact3(maxrec)
      common /c22/a,/cmf/mf,/crecord/nrecord,/cnc/nc,/csize/m
     2 ,/cic1/lact3

      do i=1,m
       grad(i)=0.0d+00
      end do     

      do i=1,nrecord
       j1=lact3(i)
       do j=1,mf
        k1=(j1-1)*mf+j
        grad(k1)=grad(k1)+2.0d+00*(x(k1)-a(i,j))
       end do
      end do
      return
      end

c=====================================================================
      subroutine dgm(x,f2)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000, maxdg=1000,maxit=100000)
      double precision x(maxvar),x1(maxvar),g(maxvar),v(maxvar)
     1 ,w(maxdg,maxvar),prod(maxdg,maxdg),z(maxdg),fvalues(maxit)
      INTEGER ij(maxdg)
      common /csize/m,/cij/ij,jvertex,/cz/z,/ckmin/kmin,/cnc/nc,/cns/ns
c====================================================================
      dist1=1.0d-07
      step0=-5.0d-02
      div=1.0d-01
      eps0=1.0d-07
      slinit=1.0d+00
      slmin=1.0d-05*slinit
      sdif=1.0d-05
      mturn=4
      maxiter=5000
      niter=0
      nbundle=min(m+3,40)
c====================================================================
      sl=slinit/div
      call fv(x,f2)
  1   sl=div*sl
      IF(sl.lt.slmin) return
      do i=1,m
       g(i)=1.0d+00/dsqrt(DBLE(m))
      end do
      nnew=0
c================================================================
   2  niter=niter+1
        IF(niter.gt.maxiter) RETURN
        nnew=nnew+1
        f1=f2
        fvalues(niter)=f1
c---------------------------------------------------------------
        if (nnew.gt.mturn) then
         mturn2=niter-mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.d+00)
         IF(ratio1.LT.sdif) GO TO 1
        end if
        if (nnew.GE.(2*mturn)) then
         mturn2=niter-2*mturn+1
         ratio1=(fvalues(mturn2)-f1)/(dabs(f1)+1.d+00)
         IF(ratio1.LT.(1.d-01*sdif)) GO TO 1
        end if
c--------------------------------------------------------------
        do ndg=1,nbundle
            call dgrad(x,sl,g,v)
            dotprod=0.d+00
            do i=1,m
             dotprod=dotprod+v(i)*v(i)
            end do
            r=dsqrt(dotprod)
            IF(r.lt.eps0) GO TO 1
            IF(ndg.eq.1) then
                         rmean=r
                         kmin=1
                         rmin=r
            END if
            IF(ndg.gt.1) then
                         rmin=dmin1(rmin,r)
                         IF(r.eq.rmin) kmin=ndg
                         rmean=((ndg-1)*rmean+r)/ndg
            END if
            toler=dmax1(eps0,dist1*rmean)
            do i=1,ndg-1
             prod(ndg,i)=0.d+00
             do j=1,m
              prod(ndg,i)=prod(ndg,i)+w(i,j)*v(j)
             end do
             prod(i,ndg)=prod(ndg,i)
            end do
            prod(ndg,ndg)=dotprod
c====================================================================
            do i=1,m
             w(ndg,i)=v(i)
            end do
            call wolfe(ndg,prod)
c================================
            do i=1,m
             v(i)=0.d+00
             do j=1,jvertex
              v(i)=v(i)+w(ij(j),i)*z(j)
             END do
            END do
c================================
            r=0.d+00
            do i=1,m
             r=r+v(i)*v(i)
            end do
            r=dsqrt(r)
            if(r.lt.toler) GO TO 1
c===========================================================
             do i=1,m
              g(i)=-v(i)/r
              x1(i)=x(i)+sl*g(i)
             end do
c===========================================================
             call fv(x1,f4)
             f3=(f4-f1)/sl
             decreas=step0*r
             if(f3.lt.decreas) then
                        call armijo(x,g,f1,f5,f4,sl,step,r)
                        f2=f5
                        do i=1,m
                         x(i)=x(i)+step*g(i)
                        end do
c                        print 3, niter,ndg,f2,sl,r
                        IF(ndg.le.2) sl=1.1d+00*sl
                        GO TO 2
c 3                      format(I6,I4,3f16.8)
             end if
         END do
c=====================================================
      go to 1
      return
      end

c==============================================================
c  Subroutines Wolfe and Equations solves quadratic
c  programming problem, to find
c  descent direction, Step 3, Algorithm 2.
c===============================================================

      subroutine wolfe(ndg,prod)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000, maxdg=1000)
      common /csize/m,/w01/a,/cij/ij,jvertex,/cz/z,/ckmin/kmin
      INTEGER ij(maxdg)
      double precision z(maxdg),z1(maxdg),a(maxdg,maxdg)
     1 ,prod(maxdg,maxdg)
      j9=0
      jmax=500*ndg
      jvertex=1
      ij(1)=kmin
      z(1)=1.d+00
c=======================================
c  To calculate X
c=======================================
 1    r=0.0d+00
      do i=1,jvertex
       do j=1,jvertex
        r=r+z(i)*z(j)*prod(ij(i),ij(j))
       end do
      end do
      IF(ndg.eq.1) RETURN
c========================================
c  To calculate <X,P_J> and J
c========================================
      t0=1.0d+12
      do i=1,ndg
        t1=0.0d+00
        do j=1,jvertex
          t1=t1+z(j)*prod(ij(j),i)
        end do
        if(t1.lt.t0) then
                     t0=t1
                     kmax=i
        end if
      end do
c========================================
c  First stopping criterion
c========================================
      rm=prod(kmax,kmax)
      do j=1,jvertex
       rm=dmax1(rm,prod(ij(j),ij(j)))
      end do
      r2=r-1.d-12*rm
      if(t0.gt.r2) RETURN
c========================================
c  Second stopping criterion
c========================================
      do i=1,jvertex
       if(kmax.eq.ij(i)) RETURN
      end do
c========================================
c Step 1(e) from Wolfe's algorithm
c========================================
      jvertex=jvertex+1
      ij(jvertex)=kmax
      z(jvertex)=0.0d+00
c========================================
 2    do i=1,jvertex
       do j=1,jvertex
        a(i,j)=1.0d+00+prod(ij(i),ij(j))
       end do
      end do
      j9=j9+1
      if(j9.gt.jmax) RETURN
      call equations(jvertex,z1)
      do i=1,jvertex
       if(z1(i).le.1.0d-10) go to 3
      end do
      do i=1,jvertex
       z(i)=z1(i)
      end do
      go to 1
  3   teta=1.0d+00
      do i=1,jvertex
       z5=z(i)-z1(i)
       if(z5.gt.1.0d-10) teta=dmin1(teta,z(i)/z5)
      end do
      do i=1,jvertex
       z(i)=(1.0d+00-teta)*z(i)+teta*z1(i)
       if(z(i).le.1.0d-10) then
                          z(i)=0.0d+00
                          kzero=i
       end if
      end do
      j2=0
      do i=1,jvertex
       IF(i.ne.kzero) then
                     j2=j2+1
                     ij(j2)=ij(i)
                     z(j2)=z(i)
       END if
      end do
      jvertex=j2
      go to 2
      return
      end

      subroutine equations(n,z1)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000, maxdg=1000)
      common /w01/a
      double precision a(maxdg,maxdg),z1(maxdg),b(maxdg,maxdg)
      do i=1,n
       do j=1,n
        b(i,j)=a(i,j)
       end do
       b(i,n+1)=1.0d+00
      end do
      do i=1,n
       r=b(i,i)
       do j=i,n+1
        b(i,j)=b(i,j)/r
       end do
       do j=i+1,n
        do k=i+1,n+1
         b(j,k)=b(j,k)-b(i,k)*b(j,i)
        end do
       end do
      end do
      z1(n)=b(n,n+1)
      do i=1,n-1
        k=n-i
        z1(k)=b(k,n+1)
        do j=k+1,n
         z1(k)=z1(k)-b(k,j)*z1(j)
        END do
      end do
      z2=0.d+00
      do i=1,n
       z2=z2+z1(i)
      end do
      do i=1,n
       z1(i)=z1(i)/z2
      end do
      return
      end

c=====================================================================
c Subroutine dgrad calculates discrete gradients
c=====================================================================
      subroutine dgrad(x,sl,g,dg)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000, maxdg=1000)
      double precision x1(maxvar),g(maxvar),x(maxvar),dg(maxvar)
      common /csize/m,/cns/ns
       do k=1,m
        x1(k)=x(k)+sl*g(k)
       end do
       if(ns.eq.1) call auxgrad(x1,dg)
       if(ns.eq.2) call funcgrad(x1,dg)
      return
      end

c===========================================================
c Line search (Armijo-type), Step 5 Algorithm 2.
c===========================================================
      subroutine armijo(x,g,f1,f5,f4,sl,step,r)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000, maxdg=1000)
      common /csize/m
      double precision x(maxvar),g(maxvar),x1(maxvar)
      step=sl
      f5=f4
  1   step=2.0d+00*step
      do i=1,m
       x1(i)=x(i)+step*g(i)
      end do
      call fv(x1,f6)
      f3=f6-f1+1.0d-02*step*r
      IF(f3.gt.0.0d+00) then
       step=step/2.0d+00
       return
      END IF
      f5=f6
      GO TO 1
      return
      end

      subroutine fv(x,f)
      implicit double precision (a-h,o-z)
      PARAMETER(maxvar=7000)
      double precision x(maxvar)
      common /cns/ns
c================================================================
      if(ns.eq.1) call auxfunc(f,x)
      if(ns.eq.2) call func(f,x)
c=================================================================
      return
      end

c==================================================================
      subroutine finresult(x)
      PARAMETER(maxvar=7000, maxrec=100000, maxclust=200, maxnft=700
     1 ,maxclass=100,maxcannot=10000, maxmust=10000)
      implicit double precision(a-h,o-z)
      double precision x(maxvar),a(maxrec,maxnft),plabel(maxclass)
     1 ,fk(maxclust),sa(maxclust),sil(maxrec),dminim(maxrec)
      integer nob(maxclass),nob1(maxclass),nel(maxclust),list1(maxrec)
     1 ,nk(maxclust,maxrec),lcan(maxcannot,2),lmust(maxmust,2)
      common /c22/a,/cnft/nft,/cmf/mf,/cnclass/nclass,/cnc/nc
     1 ,/crecord/nrecord,/clist1/list1,/cnob/plabel,nob,/cnpnout/npurity
     2 ,/cncan/ncan,lcan,/cnmust/nmust,lmust,/adminim/dminim
c========================================================================
      do i=1,nc
       nel(i)=0
      end do
      
      do j=1,nrecord
       i=list1(j)
       nel(i)=nel(i)+1
       nk(i,nel(i))=j
      end do

      f=0.0d+00
      do k=1,nrecord       
       j=list1(k)
       f1=0.0d+00
       do k1=1,mf
        f1=f1+(a(k,k1)-x(k1+(j-1)*mf))**2
       END do
       f=f+f1
      end do
c====================================================================
c Calculation of violations
c====================================================================
      mcval=0
      mmval=0
      do i=1,ncan
       i1=lcan(i,1)
       i2=lcan(i,2)
       i5=list1(i1)
       i6=list1(i2)
       if(i5.ne.i6) mcval=mcval+1
      end do
      do i=1,nmust
       i1=lmust(i,1)
       i2=lmust(i,2)
       i5=list1(i1)
       i6=list1(i2)
       if(i5.eq.i6) mmval=mmval+1
      end do
       mt=ncan+nmust
       mt1=mcval+mmval
       if(mt.gt.0) viol=1.0d+02*dble(mt1)/dble(mt)
       if(mt.eq.0) viol=0.0d+00
       if(ncan.gt.0) vcan=1.0d+02*dble(mcval)/dble(ncan)
       if(ncan.eq.0) vcan=0.0d+00
       if(nmust.gt.0) vmust=1.0d+02*dble(mmval)/dble(nmust)
       if(nmust.eq.0) vmust=0.0d+00
c====================================================================
c Purity
c====================================================================
      if (npurity.eq.1) then
       n3=0
       do i=1,nc
        do j=1,nclass
         nob1(j)=0
        end do
        do j=1,nel(i)
         j1=nk(i,j)
         do k=1,nclass
          c1=plabel(k)
          IF(a(j1,nft).EQ.c1) nob1(k)=nob1(k)+1
         end do
        end do
        n2=nob1(1)
        do k=2,nclass
         if (n2.lt.nob1(k)) then
          n2=nob1(k)
         end if
        end do
        n3=n3+n2
       end do
       purity=1.0d+02*dble(n3)/dble(nrecord)
      end if
c===================================================================
c Calculation of Davies-Bouldin (DB) validity index
c=====================================================
      do i=1,nc
       fk(i)=0.0d+00
      end do

      fdb=0.0d+00
      do i=1,nrecord
       k=list1(i)
       fk(k)=fk(k)+dminim(i)
      end do
      do i=1,nc
       IF(nel(i).gt.0) fk(i)=fk(i)/DBLE(nel(i))
      end do
      do k=1,nc
       fm=0.0d+00
       do i=1,nc
        if (i.ne.k) then
         fk2=fk(i)+fk(k)
         f12=0.0d+00
         do j=1,mf
          f12=f12+(x(j+(i-1)*mf)-x(j+(k-1)*mf))**2
         end do
         f2=fk2/f12
         fm=dmax1(fm,f2)
        end if
       end do
       fdb=fdb+fm
      end do
      db=fdb/DBLE(nc)
c============================================================
c Calculation of Dunn (Dunn) validity index
c============================================================
      if(nc.eq.1) dunn=1.0d+00
      if(nc.gt.1) then
       dn=1.0d+26
       do i=1,nc
        dn2=1.0d+26
        do j=i+1,nc
         dm=1.0d+26
         do j1=1,nel(i)
          k1=nk(i,j1)
          do j2=1,nel(j)
           k2=nk(j,j2)
           dm1=0.0d+00
           do j3=1,mf
            dm1=dm1+(a(k1,j3)-a(k2,j3))**2
           end do
           dm=dmin1(dm,dm1)
          end do
         end do
         dn2=dmin1(dn2,dm)
        end do
        dn=dmin1(dn,dn2)
       end do

       dn3=0.0d+00
       do i=1,nc
        diam=0.0d+00
        do j=1,nel(i)
         j1=nk(i,j)
         do j2=j+1,nel(i)
          j3=nk(i,j2)
          d4=0.0d+00
          do j4=1,mf
           d4=d4+(a(j1,j4)-a(j3,j4))**2
          end do         
          diam=dmax1(diam,d4)
         end do
        end do
        dn3=dmax1(dn3,diam)
       end do
       dunn=dn/dn3
      end if
c============================================================
c Calculation of silhouette
c============================================================
      ns=0
      do i=1,nrecord
       k1=list1(i)
       do k=1,nc
        sa(k)=0.0d+00
       end do
       do j=1,nrecord
        if(j.ne.i) then
         d1=0.0d+00
         do j2=1,mf
          d1=d1+(a(i,j2)-a(j,j2))**2
         end do
         k2=list1(j)
         sa(k2)=sa(k2)+d1
        end if
       end do
       da1=sa(k1)/dble(nel(k1))
       da2=1.0d+26
       do k=1,nc
        if(k.ne.k1) then
         da3=sa(k)/dble(nel(k))
         da2=dmin1(da2,da3)
        end if
       end do 
       sil(i)=(da2-da1)/dmax1(da1,da2)
       if(sil(i).gt.0.0d+00) ns=ns+1
      end do
c====================================================      
c Normalized mutual information
c====================================================
      if (npurity.eq.1) then
       hc=0.0d+00
       do i=1,nclass
        hc1=dble(nob(i))/dble(nrecord)
        hc=hc-hc1*dlog(hc1)
       end do
       ha=0.0d+00
       do i=1,nc
        ha1=dble(nel(i))/dble(nrecord)
        ha=ha-ha1*dlog(ha1)
       end do 
       hac=5.0d-01*(ha+hc)
       
       ti=0.0d+00
       do i=1,nc
        do j=1,nclass
         nob1(j)=0
        end do
        do j=1,nel(i)
         j1=nk(i,j)
         do k=1,nclass
          IF(a(j1,mf+1).EQ.plabel(k)) nob1(k)=nob1(k)+1
         end do
        end do
        
        do k=1,nclass
         if(nob1(k).gt.0) then
          ti1=dble(nob1(k))/dble(nrecord)
          ti2=dble(nrecord*nob1(k))/dble(nel(i)*nob(k))
c          ti2=dble(nob1(k))/dble(nel(i))
          ti=ti+ti1*dlog(ti2)
         end if 
        end do
       end do
       ti=ti/hac 
      end if
c======================================================================
      if(npurity.eq.2) purity=0.0d+00
      write(39,603) nc,f,db,purity,viol,vcan,vmust,ns,ti
 603  format(i5,f24.4,f9.4,f9.3,3f9.3,i10,f9.5)
      return
      end 
