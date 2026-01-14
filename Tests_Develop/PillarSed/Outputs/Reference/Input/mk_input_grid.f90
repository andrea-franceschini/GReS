program mk_input
! Built the input file: mesh and boundary conditions.
!
! nLayers  : number of initial grid layers                - xdim/ydim  : total x and y length of the rectangular base 
! nDeltax/y: number of elemental columns on the x/y axis  - dx/dy/dz   : horizontal and dvertical mesh discretization
! nbot     : number of nodes on the bottom/on a layer     - omeName    :  filename of the sedimentation rate of one material
! omega    : sed. rate, constants in space  
!******************************************************************************************************************
!
! Modules

implicit none

!
integer                                     :: i,j,nstep,nLayers,bc_flag,nX,nY,nbot,k,info,n,nDeltax,nDeltay,nneu,nn,mat,ios
integer                                     :: mm,ll,ind,ind_2,nt_row,nTriang,m,nTriang_nn,nt,np,nt_i,i_n,creep_flag
integer                                     :: strainFlag,alinfo
integer                                     :: ind_col,tetra1,tetra2,tetra3,tetra4,ind_str,n_conn,timeDiri,nummat,t
integer                                     :: omega_flag,inflag,n_connTemp,gridint
real*8                                      :: xdim, ydim, dx,dz,dy,newX,newY,time,omega,dist,elev,ome,ii
real*8                                      :: vel,xcen,ycen,ycen_dune,xcen_dune,refTime,z1,z2,z3,hmax,distmax
real*8                                      :: e0,Cr,Cc,Kz0,Ck,sigmap,sigma0,Calpha,gs,gw,beta,kxkz,multk
integer,allocatable                         :: top(:),triang_nn(:),strtop(:),ind_triang(:)
integer,allocatable                         :: triang(:,:),nod_list(:,:),tetra(:,:)
real*8,allocatable                          :: x(:),y(:),z(:),omegaVar(:,:),dirbot(:)
character(len=1000)                         :: omeName,matName,matfile,flag_dis
character(len=60)                           :: fn1,fn2,fn3,fn4,fn7,fn11,fn12,fn18,fn19
logical                                     :: ome_flag,gridFlag
!
open(13,file='input.fnames',status='old')
read(13,'(a60)') fn1
read(13,'(a60)') fn2
read(13,'(a60)') fn3
read(13,'(a60)') fn4
read(13,'(a60)') fn11
read(13,'(a60)') fn7
read(13,'(a60)')
read(13,'(a60)') fn18
read(13,'(a60)') fn19
!!
open(1,file=fn1,status='unknown')
open(2,file=fn2,status='unknown')
open(3,file=fn3,status='unknown')
open(4,file=fn4,status='unknown')
open(11,file=fn11,status='unknown')
open(7,file=fn7,status='unknown')
open(18,file=fn18,status='unknown')
open(19,file=fn19,status='unknown')
open(17,file='bath.dat',status='unknown')
open(31,file='grid_real_coord',status='unknown')
!!
!
gridFlag = .false.
read(1,*) gridint,inflag,creep_flag
if (gridint .eq. 1) gridFlag = .true.
if (gridFlag) then
   read(13,'(a60)') fn12
   open(12,file='mesh1')
   read(1,*) nLayers
else   
   read(1,*) nLayers, xdim, ydim, nDeltax,nDeltay
end if
read(1,*) dz, bc_flag,nummat
read(1,'(a10)') matfile
read(1,*) gw, beta 
open(9,file='./'//matfile//'',status='unknown')
write(9,'(i3,1x,a60)') inflag,'                           !     1-> input from advaced config, 0-> start form init'
write(9,'(e15.4,1x,a40)') gw,'                !   gw'
write(9,'(e15.4,1x,a40)') beta,'              !   beta'
write(9,'(i3,1x,a75)') creep_flag,'                     !   1-> elasoplastic, 0-> creep'
write(9,'(i3,1x,a60)') nummat,'                           !    num material'
!
if (.not. gridFlag) then
   dy = dx
   nX = nDeltax+1
   nY = nDeltay+1
   dx = xdim/float(nDeltax)
   dy = ydim/float(nDeltay)
!
   if (dx .ne. dy) then
      write(*,*) 'discretizations x and y not equal, continue? y-yes, n-no'
      read(*,*) flag_dis
      if (flag_dis .eq. 'n') stop 
   end if 
   nbot = nX*nY
   n = nbot*nLayers*(1+nLayers)
   allocate(x(n),y(n),z(n),stat=info)
!   write(*,*) 'n ',n
   !! coordinate 2d grid
   k = 1
   do i = 1,nY
      do j = 1,nX
         x(k) = (j-1)*dx
         y(k) = (i-1)*dy
         z(k) = 0d0
         k = k+1
      end do   
   end do
   nTriang = (nX-1)*(nY-1)*2
   nt_row = nX*2-2
else   

   read(12,*) nTriang,n 
   nbot = n
   allocate(x(n),y(n),z(n),stat=info)
!
!
   if (gridFlag) then
      do j = 1,n
         read(12,*) k,x(j),y(j) 
         x(j)=x(j)!*0.001d0 
         y(j)=y(j)!*0.001d0 
         read(17,*) z1,z2,z(j)
         write(*,*) z(j)
         !z(j)=-z(j)
      end do
   !!! bathymetry aert
   !
   end if
   n = nbot*nLayers*(1+nLayers)
end if
!
!! 2d topology matrix

allocate(triang(nTriang,4))
if (.not. gridFlag) then
  do j = 1,nY-1
      mm = 0
      ll = 0
      ind = (j-1)*nt_row
      ind_2 = (j-1)*nX
      do i = 1,nt_row
         k = ind+i
      ! element 1 -sx
         if (mod(k,2) .eq. 1) then
            triang(k,2) = ind_2+i-mm;
            triang(k,3) =  ind_2+i-mm+1;
            triang(k,4) =  ind_2+i-mm+nX;
            mm = mm + 1;
         else
            triang(k,2) =  ind_2+i-ll;
            triang(k,3) =  ind_2+i+nX-ll;
            triang(k,4) =  ind_2+i+nX-ll-1;
            ll = ll + 1;
         end if
      end do 
   end do   
   n_conn = 2*1+2*2+(nX-2)*6+(nY-2)*6+(nbot-nx*2-nY*2+4)*6
   write(3,300)nTriang,n_conn
   !write(*,*)'nTriang ',nTriang,n_conn
   allocate(triang_nn(n_conn))
else
   do i = 1,nTriang
      read(12,*) triang(i,1),(triang(i,j),j=2,4)
   end do   
end if
!
allocate(top(nbot))
do i =1,nbot
   top(i)  = i
end do  
! 2D topology
allocate(ind_triang(nbot))
m=1;
if (gridFlag) then
   n_conn = 0
   n_connTemp = n*10
   allocate(triang_nn(n_connTemp))
end if    
do i=1,nbot
   ind_triang(i)=m;
   do j=1,nTriang
      do k=2,4
         if (i==triang(j,k)) then
         !   if (.not. gridFlag)  write(3,*) j
            if (gridFlag)  n_conn = n_conn+1  
              triang_nn(m)=j
              m=m+1
         endif
      end do
   end do
end do
!
if (gridFlag)  write(3,*)nTriang,n_conn
do i = 1,nTriang
   write(3,*) (triang(i,j),j=2,4)
   !write(*,*) 'triang',(triang(i,j),j=2,4)
end do
do i = 1,n_conn
   write(3,*) triang_nn(i)
end do
!allocate(triang_nn(n_conn))
!do i = 1,n_conn
!   triang_nn(i) = 
!deallocate(triang_nn)
!
do i = 1,nbot
   write(3,*) ind_triang(i)
   !write(*,*) ind_triang(i)
end do
!
allocate(nod_list(nLayers+1,nbot),strtop(nbot))
do k =1,nbot
   nod_list(1,k) = top(k)
   strtop(k) = 1
end do
!************************************************
! 3D coordinates
nt = n_conn*nLayers
allocate(tetra(nt,4))
write(2,'(3x,f5.3)') dz
write(2,300) nbot,nLayers,bc_flag
write(2,300) n,n_conn
write(31,*) nbot
!!!! scaled
do j=1,nLayers+1
   do i =1,nbot
      !write(2,100) i+(j-1)*nbot,x(i),y(i),dz*(j-1)+z(i)
      write(2,100) i+(j-1)*nbot,(x(i)-x(1))/1.0d0,(y(i)-y(1))/1.0d0,dz*(j-1)+z(i)
      write(31,100) i+(j-1)*nbot,x(i),y(i),dz*(j-1)+z(i)
    end do
end do
!
nt_i = 0
i_n=nbot
nTriang_nn = 0
do i=1,nLayers
   do ind_col=1,nbot
      i_n = i_n + 1;
      ind_str=strtop(ind_col)
      nod_list(ind_str+1,ind_col) = i_n
      top(ind_col) = i_n
      strtop(ind_col) = ind_str + 1
      newX=x(nod_list(ind_str,ind_col))
      newY=y(nod_list(ind_str,ind_col))
      !! Topologia 3D
      ! check all the edges
      if (ind_col==nbot) then
         nTriang_nn=n_conn-ind_triang(ind_col)+1;
      else
         nTriang_nn=ind_triang(ind_col+1)-ind_triang(ind_col);
      end if
      do  k=1,nTriang_nn
         nt_i = nt_i+1
         tetra(nt_i,1) = nod_list(ind_str,ind_col);
         tetra(nt_i,2) = top(triang(triang_nn(ind_triang(ind_col)+(k-1)),2));
         tetra(nt_i,3) = top(triang(triang_nn(ind_triang(ind_col)+(k-1)),3));
         tetra(nt_i,4) = top(triang(triang_nn(ind_triang(ind_col)+(k-1)),4));
         !write(*,*)  k,(tetra(nt_i,j),j=1,4),ind_col,'w',ind_triang(ind_col),'q',triang_nn(ind_triang(ind_col))
      enddo
   end do
end do
!
do i=1,nt
   write(2,200) i,(tetra(i,j),j=1,4),1
enddo
!
!
!! Dirichlet b.c.
k=1
!if (bc_flag == 0) then
!   do i=1,n
!      if coord(i,2)==coord(1,2)
!        contp(k)=i;
!        k=k+1;
!      end
!   end
!   contp=[contp top];
!   contp=unique(contp);
!   np=size(contp,2);
!else
!   contp=top;
np=nbot;
!end
!
write(4,*) np
write(4,'(6(3x,i15))') (top(i),i=1,nbot)
timeDiri = 0d0
write(4,'(i5)') timeDiri
!dirvalue=zeros(np,1);
!do i=1,np
!   dirvalue(i) = 0.0;
!write(*,*) 'test con Dir bottom = 3'
write(4,'(20(2x,f4.2))') (0.d0,i=1,nbot)

!!!write(19,*) nbot
!!!allocate(dirbot(nbot),stat=info)
!!!write(19,'(6(3x,i15))') (i,i=1,nbot)
!!!do i = 1,nbot
!!!  read(18,*) ii,ii,dirbot(i)
!!!  !write(19,'(6(3x,i15))') i
!!!end do
!!!timeDiri = 0d0
!!!write(19,'(i5)') timeDiri
!!!write(19,'(10(2x,f9.2))') (dirbot(i),i=1,nbot)

!enddo
!**********************
nneu = 0
write(11,'(i8)') nneu
!
do i =1,nbot
   write(7,*) top(i),2
end do
!
do i =1,nbot
   write(7,*) (nod_list(j,i),j=1,2)
end do
!  material 
do i=1,nummat
   read(1,*) matName
   open(10,file='./materials/'//matName//'',status='unknown')
   write(9,'(a60)') './Input/materials/'//matName//''
   read(1,*) e0,Cr,Cc,Kz0,Ck
   read(1,*) multk,sigmap,sigma0,Calpha,gs,refTime
!   write(*,*) multk,sigmap,sigma0,Calpha,gs,refTime
   write(10,500) e0,Cr,Cc,Kz0,multk,Ck,sigmap,sigma0,Calpha,gs,refTime
   write(*,*) 'mat ',nummat,'(Cc-Cr)/(1+e0)/Calpha',(Cc-Cr)/(1+e0)/Calpha
   if ((Cc-Cr)/(1+e0)/Calpha .gt. 100d0) then
      write(*,*) '(CR-RR)/Calpha > 100!! check parameters'     
      write(*,*) 'continue? y-yes, n-no'
      read(*,*) flag_dis
      if (flag_dis .eq. 'n') stop
   end if
   write(10,'(a60)') '  e0          Cr         Cc        Kz0        Kx/Kz '
   write(10,'(a60)') '  Ck          sigmap     sigma0    Calpha     gs    '
   write(10,'(a60)') '  reference time    '
   close(10)
end do
! omega
read(1,*) omega_flag
ome_flag = .false.
do i=1,nummat
   read(1,*) omeName
   write(9,'(a60)') './Input/materials/'//omeName//''
   open(20+i,file='./materials/'//omeName//'',status='unknown')
end do
!
!   omega - > 0: do nothing - 1: print constant omega - 2: lobe examples
if (omega_flag .eq. 0)  write(*,*) 'omegafile not overwrite'
if (omega_flag .eq. 1)  then
   write(*,*) 'make omegafile with constant values'  
   do i=1,nummat
!      open(10,file='./materials/'//omeName//'',status='unknown')
      read(1,*,IOSTAT=ios) nstep
      if (ios .gt. 0) then
           write(*,*) 'wrong input nstep'
           stop
      end if
      !write(*,*) nstep
      do k = 1,nstep
        read(1,*,IOSTAT=ios) time,omega
        if (ios .gt. 0) then
           write(*,*) 'wrong input for omegafile'
           stop
        end if   
        write(*,*) 'time',time,omega
        write(20+i,'(3x,f12.6)') time
        write(20+i,'(6(1x,e15.5))') (omega,j=1,nbot)
      end do
      close(20+i)
   end do
!
elseif (omega_flag .eq. 2) then
   write(*,*) 'omegafile -> delta lobe case'  
   allocate(omegaVar(nummat,nbot))
   ycen = 5d0
   vel = 1/100d0
   elev = 9d0/4d0
!   do i=1,nummat
!      open(10+i,file='./materials/'//omeName//'',status='unknown')
!   end do
   do t = 1,1000,1
!
      time = float(t)
      xcen = (time-400d0)*vel
      !write(*,*) 'xcen',xcen
      do i=1,nbot
         !nn=top(i)
         dist = sqrt((x(i)-xcen)**2d0+(y(i)-ycen)**2d0)
         write(*,*)x(i),dist,time,y(i)
         omegaVar(1,i) =  0.0
         omegaVar(2,i) =  0.0
         omegaVar(3,i) =  0.0
!
         if (dist .lt. 4d0) then
            if (dist .lt. 4.0d0 .and. dist .ge. 2.41d0) omegaVar(1,i) = 2d0*elev*vel*(abs(x(i)-xcen))/dist
            !if (y(i) = ycen) write(111,*) omegaVar(1,i)
            !write(*,*) omegaVar(1,i),time
            if (dist .lt. 2.41d0 .and. dist .gt. 1.81d0) omegaVar(2,i) = 2d0*elev*vel*(abs(x(i)-xcen))/dist
            if (dist .lt. 1.81d0 .and. dist .gt. 1.5d0) omegaVar(3,i) = 2d0*elev*vel*(abs(x(i)-xcen))/dist
            if (x(i) .lt. xcen) then
               do mat=1,3
                  omegaVar(mat,i) = 0.0
               end do   
!         else
!            omega(i)=0d0
            end if
         end if
!    omega(i)=omega(i)*2d0      
!!!   
         ome = 0d0
         if ((time .le. 600 .and. time .ge. 500.) .or. (time .le. 1000 .and. time .ge. 900.)) then
            if (time .le. 600 .and. time .ge. 500.) xcen_dune = 2d0
            if (time .le. 1000 .and. time .ge. 900.) xcen_dune = 6d0
            if (y(i) .ge. 3.5d0 .and. y(i) .le. 6.5d0)  then
               if (y(i) .ge. 4d0 .and. y(i) .le. 6d0)  then
                  ome = -32d0*(x(i)-xcen_dune)**2+4d0
                  ome = ome/100d0
        !         mat_time(i) = 3
                  if (ome .lt. 0) ome = 0
               elseif (y(i) .ge. 6d0 .and. y(i) .le. 6.5d0)  then
                  ycen_dune = 6d0
                  ome = -32d0*(x(i)-xcen_dune)**2-16d0*(y(i)-ycen_dune)**2+4d0
                  ome = ome/100d0
        !         mat_time(i) = 3
                  if (ome .lt. 0) ome = 0
               elseif (y(i) .ge. 3.5d0 .and. y(i) .le. 4d0)  then
                  ycen_dune = 4d0
                  ome = -32d0*(x(i)-xcen_dune)**2-16d0*(y(i)-ycen_dune)**2+4d0
                  ome = ome/100d0
        !         mat_time(i) = 3
                  if (ome .lt. 0) ome = 0
               end if
            end if
            if (ome .gt. 0) omegaVar(3,i) = ome + omegaVar(3,i)
         end if
!      !!!!    
      end do
! write here!!
   ! print time and omega
      do k = 1,nummat
         write(20+k,'(f7.2)') time-1d0
         write(20+k,'(6(1x,e15.5))') (omegaVar(k,j),j=1,nbot)
      end do
   end do
!   
   do i=1,nummat
      close(20+i)
   end do   
!
elseif (omega_flag .eq. 3) then
   write(*,*) 'omegafile -> meander oxbow case'
   write(*,*) x
   write(*,*) y
   allocate(omegaVar(nummat,nbot),stat=alinfo)
   if (alinfo .ne. 0) stop

      do i=1,nbot
         do j=1,nummat
            omegaVar(j,i)=0.0
         end do
      end do
   !******************************! MEANDER 
   hmax = 3.d0
   distmax = 5d0
   do t = 0,2200,1
      do i=1,nbot
         do j=1,nummat
            omegaVar(j,i)=0.0
         end do
      end do
      time = float(t)
      !write(*,*) time
      do i=1,nbot
         !nn=top(i)
         nn=i
         dist=sqrt((y(nn)-distmax)**2+x(nn)**2)   ! distance nn - (5,0)
         write(*,*) dist,x(nn),y(nn),nn
         ! formazione del meandro time < 100
         if (time .le. 100d0) then
            if (dist .lt. 5d0) ome=hmax
         !!!      if (dist .lt. distmax) omega(i)=hmax+(dist/distmax-1d0)*0.5d0 
            if (dist .ge. distmax .and. dist .lt. 5.85d0)      ome=hmax-((dist-distmax)/0.85d0*2d0)
            if (dist .ge. 5.85d0 .and. dist .lt. 7.15d0)   ome=abs(dist-6.50d0)/0.65d0*0.5d0+0.5d0
            if (dist .ge. 7.15d0 .and. dist .lt. 8d0)      ome=hmax-((8d0-dist)/0.85d0*2d0)
            if (dist .ge.  8d0)   ome=hmax
            ome = ome/100d0
            omegaVar(1,i) = ome
            write(*,*) ome,time,dist,x(nn),x(nn)
            !! infill  
         elseif (time .ge. 100d0 .and. time .lt. 2100d0) then
                 !! every 200 years the infill material changes
         ! 
            if (dist .lt. 5d0 .or. dist .ge. 8d0) then
             !     omega(i) = 0.5d0
               ome = 0.0d0
             !    omega(i) = ((5d0-dist)/10d0+0.5d0)*0.6d0 
               omegaVar(2,i) = ome
               else if (dist .ge.  5d0 .and. dist .lt. 5.85d0) then
                   ome=(dist-5d0)/0.85d0*2.0d0
   !
                  else if (dist .ge.  7.15d0 .and. dist .lt. 8.d0) then
                   ome=(8.0d0-dist)/0.85d0*2.0d0
                  else if (dist .ge.  5.85d0 .and. dist .lt. 7.15d0) then
                ome=2.0d0+(0.65d0-abs(dist-6.5d0))/0.65d0*0.3d0
   !
            else
               ome = 0d0
            end if
               ome = ome/2000d0
               !
            if (mod(int(time/200d0),2) .eq. 0) then
               omegaVar(2,i) = ome
            else
                  omegaVar(3,i) = ome
            end if
            !!!!! from 2100 to 2200 
        else
            ome = 0.5d0/100d0
            omegaVar(4,i) = ome
         end if
!
! write here!!
   ! print time and omega
      end do
      do k = 1,nummat
         write(20+k,'(f7.2)') time
         write(20+k,'(6(1x,e15.5))') (omegaVar(k,j),j=1,nbot)
      end do
!
   end do
   do i=1,nummat
      close(20+i)
   end do
!
end if        
!
!
100 format(i8,3(3x,e19.12))
200 format(6(3x,i8))
300 format(3(3x,i8))
400 format(4(3x,e19.12))
500 format(5(3x,e19.12))
end program mk_input
