program main

implicit none

integer, parameter :: natom=2,nxyz_x=3,nxyz_y=6,nxyz_z=5
integer, parameter :: nA=nxyz_x*nxyz_y*nxyz_z
integer, parameter :: nB=nxyz_x*nxyz_y*nxyz_z
integer, parameter :: N_pos=natom*nxyz_x*nxyz_y*nxyz_z
integer :: nx,ny,nz,n,m,sites
integer :: i,j,k
integer :: NL1st(N_pos,N_pos),NN1st(N_pos)
                                        ! otherwise, codes has to be changed

double precision :: a,b,c,latt(3),dist_tmp(3),length
double precision :: pos(N_pos,3),cellpos(natom,3),dist(N_pos,N_pos,3)
double precision :: rcutAA1st,rcutBB1st,rcutAB1st,rcutBA1st
double precision :: rcut_lower1st
!!!!!!!Cubic perovskites with a,b,c and atomic positions!!!!!
a=3.8
b=3.8
c=3.8
cellpos(1,1)=0.0
cellpos(1,2)=0.0
cellpos(1,3)=0.0
cellpos(2,1)=0.5
cellpos(2,2)=0.5
cellpos(2,3)=0.5
rcutAA1st=a*1.01   ! should be slightly larger than a
rcutBB1st=a*1.01
rcutAB1st=a*1.01/sqrt(2.0)
rcutBA1st=a*1.01/sqrt(2.0)
rcut_lower1st=0.01
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!Generating the atomic positions in a simulated box and output the
!!!!!!initial-poscars for check!!!!!!
call init_pos(natom,N_pos,nxyz_x,nxyz_y,nxyz_z,a,b,c,cellpos,pos,latt)
call out_init_poscar(natom,N_pos,nxyz_x,nxyz_y,nxyz_z,pos,latt)


!!!!!!!computing the distance between two atoms, and save it in an array
open(11,file='_distance')
write(11,*)"initial distance of i and j atoms"


do i=1,N_pos
  do j=1,N_pos
       call distance(i,j,nxyz_x*a,nxyz_y*b,nxyz_z*c,N_pos,pos,dist_tmp)
    do k=1,3
       dist(i,j,k)=dist_tmp(k)
    end do
       length=sqrt(dist(i,j,1)*dist(i,j,1)+dist(i,j,2)*dist(i,j,2)+dist(i,j,3)*dist(i,j,3))
       write (11,"(2I4.0,4f12.6)") i,j,dist(i,j,1),dist(i,j,2),dist(i,j,3),length
  enddo
end do

close(11)

!!!!!!find the 1st nearest nrighbor of each atom
NL1st(:,:)=0
NN1st(:)=0 ! 12 for A-B, 21 for B-A, 11 for A-A, 22 for B-B, etc


open(11,file='_1stNeighborAA')
sites=11
call neighbor(sites,N_pos,nA,nB,rcut_lower1st,rcutAA1st,dist,NN1st,NL1st)
do i=1,N_pos
   do j=1,N_pos
      if (NL1st(i,j) .ne. 0) then
        write(11,'(3I6.0)') i,NN1st(i),NL1st(i,j)
      end if
   end do
end do 
close(11)

open(11,file='_1stNeighborBB')
sites=22
call neighbor(sites,N_pos,nA,nB,rcut_lower1st,rcutAA1st,dist,NN1st,NL1st)
do i=1,N_pos
   do j=1,N_pos
      if (NL1st(i,j) .ne. 0) then
        write(11,'(3I6.0)') i,NN1st(i),NL1st(i,j)
      end if
   end do
end do
close(11)

open(11,file='_1stNeighborAB')
sites=12
call neighbor(sites,N_pos,nA,nB,rcut_lower1st,rcutAA1st,dist,NN1st,NL1st)
do i=1,N_pos
   do j=1,N_pos
      if (NL1st(i,j) .ne. 0) then
        write(11,'(3I6.0)') i,NN1st(i),NL1st(i,j)
      end if
   end do
end do
close(11)

open(11,file='_1stNeighborBA')
sites=21
call neighbor(sites,N_pos,nA,nB,rcut_lower1st,rcutAA1st,dist,NN1st,NL1st)
do i=1,N_pos
   do j=1,N_pos
      if (NL1st(i,j) .ne. 0) then
        write(11,'(3I6.0)') i,NN1st(i),NL1st(i,j)
      end if
   end do
end do
close(11)



end 




!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!
! cccccccc subroutine 1: generating the initial positions of atoms
subroutine init_pos(natom,N_pos,nxyz_x,nxyz_y,nxyz_z,a,b,c,cellpos,pos,latt)
  integer natom,nxyz_x,nxyz_y,nxyz_z
  integer nx,ny,nz,m,n  
  double precision a,b,c,cellpos(natom,3),pos(N_pos,3),latt(3)
  
  n=1
  pos(:,:)=0.0d0
  latt(1)=a*nxyz_x
  latt(2)=b*nxyz_y
  latt(3)=c*nxyz_z
   do m=1,natom
     do nx=0,nxyz_x-1
       do ny=0,nxyz_y-1
         do nz=0,nxyz_z-1
           pos(n,1)=a*nx+a*cellpos(m,1)
           pos(n,2)=b*ny+b*cellpos(m,2)
           pos(n,3)=b*nz+c*cellpos(m,3)
           n=n+1
         end do
       end do
     end do
   end do

end subroutine

!ccccccc subroutine 2: output the generated initial positions
subroutine out_init_poscar(natom,N_pos,nxyz_x,nxyz_y,nxyz_z,pos,latt)
  integer natom,N_pos,nxyz_x,nxyz_y,nxyz_z
  integer nx,ny,nz,m,n,nA,nB
  double precision pos(N_pos,3),latt(3)

  nA=1*nxyz_x*nxyz_y*nxyz_z
  nB=1*nxyz_x*nxyz_y*nxyz_z
  nO=3*nxyz_x*nxyz_y*nxyz_z
  n=1

     open(11,file='_poscar.vasp') 
     write(11,*)"Perovskites"
     write(11,'(f12.8)') 1.0
     write(11,'(3f12.8)') latt(1),0.000,0.000
     write(11,'(3f12.8)') 0.000,latt(2),0.000
     write(11,'(3f12.8)') 0.000,0.000,latt(3)
     write(11,*)"Ba  Ti"
     write(11,"(3I4.2)") nA,nB
     write(11,*)"Cartesian"

     do m=1, natom
       do nx=0, nxyz_x-1
         do ny=0, nxyz_y-1
           do nz=0, nxyz_z-1
           write(11,'(3f12.8)') pos(n,1),pos(n,2),pos(n,3)
           n=n+1
           end do
         end do
       end do
     end do

     close(11)

end subroutine

!cccccccc subroutine 3: distance between atoms i and j, periodic boundary
!condictions is considered
subroutine distance(i,j,La,Lb,Lc,N_pos,pos,d_ij)
  integer i,j,N_pos,n
  double precision La,Lb,Lc,pos(N_pos,3)
  double precision d_ij(3),L(3)

 L(1)=La
 L(2)=Lb
 L(3)=Lc

  do n=1,3
    d_ij(n)=pos(i,n)-pos(j,n)
    if (d_ij(n)<(-1.0*L(n)/2.0)) then
        d_ij(n)=d_ij(n)+L(n)
    else if (d_ij(n)>(L(n)/2.0)) then
        d_ij(n)=d_ij(n)-L(n)
    else
        d_ij(n)=d_ij(n)
    end if
  end do

end subroutine
  
!ccccccc subroutine 4: to find the nearest neighbours of each atom
subroutine neighbor(sites,N_pos,nA,nB,rcut1,rcut2,dist,NN,NL)
  integer N_pos,NN(N_pos),NL(N_pos,N_pos)
  integer sites,i,j,nA,nB
  double precision rcut1,rcut2,dist(N_pos,N_pos,3),d_square
  double precision rcut_square1,rcut_square2

  NN(:)=0
  NL(:,:)=0
  rcut_square1=rcut1*rcut1
  rcut_square2=rcut2*rcut2
 
  if (sites .eq. 11) then ! 11 for A-A, 12 for A-B, 21 for B-A, 22 for B-B,etc
    do i=1,nA
      do j=1,nA
          d_square=dist(i,j,1)*dist(i,j,1)+dist(i,j,2)*dist(i,j,2)+dist(i,j,3)*dist(i,j,3)
            if ((d_square .lt. rcut_square2) .and. (d_square .gt. rcut_square1)) then
                NN(i)=NN(i)+1
                NL(i,NN(i))=j
            end if
      end do
    end do

  else if (sites .eq. 22) then
    do i=nA+1,nA+nB
      do j=nA+1,nA+nB
          d_square=dist(i,j,1)*dist(i,j,1)+dist(i,j,2)*dist(i,j,2)+dist(i,j,3)*dist(i,j,3)
            if ((d_square .lt. rcut_square2) .and. (d_square .gt. rcut_square1)) then
                NN(i)=NN(i)+1
                NL(i,NN(i))=j
            end if
      end do
    end do

  else if (sites .eq. 12) then
    do i=1,nA
      do j=nA+1,nA+nB
          d_square=dist(i,j,1)*dist(i,j,1)+dist(i,j,2)*dist(i,j,2)+dist(i,j,3)*dist(i,j,3)
            if ((d_square .lt. rcut_square2) .and. (d_square .gt. rcut_square1)) then
                NN(i)=NN(i)+1
                NL(i,NN(i))=j
            end if
      end do
    end do

  else if (sites .eq. 21) then
    do i=nA+1,nA+nB
      do j=1,nA
          d_square=dist(i,j,1)*dist(i,j,1)+dist(i,j,2)*dist(i,j,2)+dist(i,j,3)*dist(i,j,3)
            if ((d_square .lt. rcut_square2) .and. (d_square .gt. rcut_square1)) then
                NN(i)=NN(i)+1
                NL(i,NN(i))=j
            end if
      end do
    end do

  else
     write(*,*)"Not considered yet!!!"
  end if

end subroutine
  




