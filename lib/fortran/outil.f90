module outil

  use declaration

contains

  !************************************************************************
  ! 									*
  !FROM MODULE  AGRIF_MODUTIL FROM NESTING_TOOLS F.LEMARIE
  !MODIFIED by Matthieu Caillaud 18/11/08	
  !                       						*
  !									*
  ! module containing subroutine used for : 				*
  !   - unrolling 2D arrays to 1D arrays (required for SCRIP package use)	*
  !   - convert 1D arrays to 2D arrays (required for SCRIP package use)	*
  !   - remapping process (use SCRIP remapping matrix)			*
  !									*
  !************************************************************************
  ! 
  
  !************************************************************************
  !   SUBROUTINE 1Dto2D
  !************************************************************************
  
  subroutine tab1Dto2D(tab1D,tab2D,nx,ny)
    
    implicit none   

  !***********************************************
  
    real(kind=rlg),dimension(:) :: tab1D
    real(kind=rlg),dimension(:,:) :: tab2D
    integer :: xpos,ypos
    integer :: nx,ny   
    integer :: i      

  !************************************************      
    xpos=0
    ypos=1
   
    do i=1,nx*ny
       xpos=xpos+1
       if(xpos.gt.nx) then
          xpos=1
          ypos=ypos+1
       endif
       tab2D(ypos,xpos)=tab1D(i)
    End Do
    
  end subroutine tab1Dto2D

  !************************************************************************
  !   END SUBROUTINE 1Dto2D
  !************************************************************************

  !************************************************************************
  !   SUBROUTINE tab2Dto1D
  !************************************************************************
  
  subroutine tab2Dto1D(tab2D,tab1D)

    implicit none 
    
    real(kind=rlg),dimension(:,:) :: tab2D
    real(kind=rlg),dimension(:) :: tab1D
    
    integer :: xpos,ypos
    integer :: nx,ny
    integer :: i
    
    nx = size(tab2D,2)
    ny = size(tab2D,1)
   
    xpos = 0
    ypos = 1
    Do i = 1,nx*ny
       xpos = xpos + 1
       if(xpos.gt.nx) then
          xpos = 1
          ypos = ypos + 1
       end if
       tab1D(i) = tab2D(ypos,xpos)
    End Do

  end subroutine tab2Dto1D

  !************************************************************************
  !   END SUBROUTINE tab2Dto1D
  !************************************************************************

  
  !************************************************************************
  !   SUBROUTINE tab2Dto1D logical
  !************************************************************************


  subroutine logtab2Dto1D(tab2D,tab1D)

    implicit none 
    
  !*********************************************       
  
    logical,dimension(:,:) :: tab2D
    logical,dimension(:) :: tab1D
    integer :: xpos,ypos
    integer :: nx,ny
    integer :: i

  !**********************************************
  
    nx = size(tab2D,2)
    ny = size(tab2D,1)
        
    xpos = 0
    ypos = 1
    Do i = 1,nx*ny
       xpos = xpos + 1
       if(xpos.gt.nx) then
          xpos = 1
          ypos = ypos + 1
       end if
       tab1D(i) = tab2D(ypos,xpos)
    End Do
 
  end subroutine logtab2Dto1D

  !************************************************************************
  !   END SUBROUTINE tab2Dto1D logical
  !************************************************************************ 


  !************************************************************************
  !   SUBROUTINE 1Dto2D logical
  !************************************************************************

  subroutine logtab1Dto2D(tab1D,tab2D,nx,ny)

    implicit none   

  !****************************************
  
    logical,dimension(:) :: tab1D
    logical,dimension(:,:) :: tab2D
    integer :: xpos,ypos
    integer :: nx,ny   
    integer :: i      

  !****************************************      
  
    xpos=0
    ypos=1
    
    do i=1,nx*ny
       xpos=xpos+1
       if(xpos.gt.nx) then
          xpos=1
          ypos=ypos+1
       endif
       tab2D(ypos,xpos)=tab1D(i)
    End Do
    
  end subroutine logtab1Dto2D

  !************************************************************************
  !   END SUBROUTINE 1Dto2D logical
  !************************************************************************
  
  !**************************************************************
  !   subroutine make_remap
  !**************************************************************            

  subroutine make_remap(tabin,tabout,nxfin,nyfin,matrix,src_add,dst_add) 
    
    implicit none

  !*************************************************************      
  
    Real(kind=rlg), dimension(:,:)  :: tabin 
    Real(kind=rlg), dimension(:,:) :: tabout 
    Real(kind=rlg), pointer, dimension(:,:) :: tabtemp    
    Integer,dimension(:) :: src_add,dst_add
    Integer :: nxfin,nyfin      
    Real(kind=rlg), pointer, dimension(:) :: var1D,var_interp1D
    Real(kind=rlg),dimension(:,:) :: matrix 
    Integer :: num_links,i

  !*************************************************************
  
    allocate(var1D(size(tabin,1)*size(tabin,2)))     
    call tab2Dto1D(tabin,var1D)

    allocate(var_interp1D(nxfin*nyfin))
    var_interp1D = 0.0
    num_links = size(dst_add)
    
    Do i = 1,num_links
       var_interp1D(dst_add(i)) = var_interp1D(dst_add(i)) &
            + matrix(1,i)*var1D(src_add(i))
    End do
    
    allocate(tabtemp(size(tabout,1),size(tabout,2)))

    Call tab1Dto2D(var_interp1D,tabtemp,nyfin,nxfin)

    tabout = tabtemp

    deallocate(var_interp1D,var1D,tabtemp)  

  end subroutine make_remap

  !**************************************************************
  !   end subroutine make_remap
  !**************************************************************          

  !**************************************************************
  !   subroutine make_conserv_remap
  !**************************************************************            

  subroutine make_conserv_remap(tabin,tabout,nxfin,nyfin,matrix,src_add,dst_add)
    
    implicit none
    
    Real(kind=rlg), dimension(:,:)  :: tabin
    Real(kind=rlg), dimension(:,:) :: tabout
    Real(kind=rlg), pointer, dimension(:,:) :: tabtemp
    Integer,dimension(:) :: src_add,dst_add
    Integer :: nxfin,nyfin
    Real(kind=rlg), pointer, dimension(:) :: var1D,var_interp1D
    Real(kind=rlg),dimension(:,:) :: matrix
    Integer :: num_links,i

  !***************************************************************
  
    allocate(var1D(size(tabin,1)*size(tabin,2)))
    call tab2Dto1D(tabin,var1D)
   
    allocate(var_interp1D(nxfin*nyfin))
    var_interp1D = 0.0
    num_links = size(dst_add)
    
    Do i = 1,num_links
       var_interp1D(dst_add(i)) = var_interp1D(dst_add(i)) &
            + matrix(1,i)*var1D(src_add(i)) 
       !  +matrix(2,i)*var1D(src_add(i))  &
       ! +matrix(3,i)*var1D(src_add(i))

    End do

    allocate(tabtemp(size(tabout,1),size(tabout,2)))

    Call tab1Dto2D(var_interp1D,tabtemp,nyfin,nxfin)

    tabout = tabtemp

    deallocate(var_interp1D,var1D,tabtemp)
    
  end subroutine make_conserv_remap
           

  !**************************************************************
  !   end subroutine make_conserv_remap
  !**************************************************************



  !**************************************************************
  !   subroutine make_dstwgt_remap
  !**************************************************************            

  subroutine make_dstwgt_remap(tabin,tabout,nxfin,nyfin,matrix,src_add,dst_add)
    
    implicit none

  !**************************************************************      
  
    Real(kind=rlg), dimension(:,:)  :: tabin
    Real(kind=rlg), dimension(:,:) :: tabout
    Real(kind=rlg), pointer, dimension(:,:) :: tabtemp
    Integer,dimension(:) :: src_add,dst_add
    Integer :: nxfin,nyfin
    Real(kind=rlg), pointer, dimension(:) :: var1D,var_interp1D
    Real(kind=rlg),dimension(:,:) :: matrix
    Integer :: num_links,i

  !**************************************************************** 
  
    allocate(var1D(size(tabin,1)*size(tabin,2)))
    call tab2Dto1D(tabin,var1D)
   
    allocate(var_interp1D(nxfin*nyfin))
    var_interp1D = 0.0
    num_links = size(dst_add)

    Do i = 1,num_links
       var_interp1D(dst_add(i)) = var_interp1D(dst_add(i)) &
            + matrix(1,i)*var1D(src_add(i))
    End do

    allocate(tabtemp(size(tabout,1),size(tabout,2)))
    
    Call tab1Dto2D(var_interp1D,tabtemp,nyfin,nxfin)

    tabout = tabtemp
    
    deallocate(var_interp1D,var1D,tabtemp)
    
  end subroutine make_dstwgt_remap
  
  !**************************************************************
  !   end subroutine make_dstwgt_remap
  !**************************************************************

  
  !**************************************************************
  !   subroutine make_bicubic_remap
  !**************************************************************            

  subroutine make_bicubic_remap(tabin,masksrc,tabout,nxfin,nyfin,matrix,src_add,dst_add) 
    
    implicit none

  !********************************************************************************************      
  
    Real(kind=rlg), dimension(:,:)  :: tabin
    Logical, dimension(:,:)  :: masksrc
    Logical, pointer, dimension(:)  :: grid1_mask
    Real(kind=rlg), dimension(:,:) :: tabout     
    Integer,dimension(:) :: src_add,dst_add
    Integer :: nxfin,nyfin      
    Real(kind=rlg), pointer, dimension(:) :: var1D,var_interp1D,gradi,gradj,gradij,deriv1,deriv2
    Real(kind=rlg),dimension(:,:) :: matrix 
    Integer :: num_links,i,j,nx,ny,n,ip1,im1,jp1,jm1
    Integer :: in,is,ie,iw,ine,inw,ise,isw
    Real(kind=rlg) :: delew,delns

  !**********************************************************************************************
  
    nx = size(tabin,1)
    ny = size(tabin,2)
    Allocate(gradi(nx*ny),gradj(nx*ny),gradij(nx*ny),deriv1(nx*ny),deriv2(nx*ny))
    Allocate(var1D(nx*ny),grid1_mask(nx*ny))
    
    Call tab2Dto1D(tabin,var1D)
    Call logtab2Dto1D(masksrc,grid1_mask)     
    
    gradi  = 0.0
    gradj  = 0.0
    gradij = 0.0

    Do n = 1,nx*ny
       
       IF( grid1_mask(n) ) then                        
          
          delew = 0.5      
          delns = 0.5
          
          j = (n-1)/ny + 1
          i = n - (j-1)*ny
          
          ip1 = i+1
          im1 = i-1
          jp1 = j+1
          jm1 = j-1      
          
          If (ip1 > ny) ip1 = ip1 - ny
          
          If (im1 < 1 ) im1 = ny
          
          If (jp1 > nx) then
             jp1 = j
             delns = 1.
          Endif
             
          If (jm1 < 1 ) then
             jm1 = j
             delns = 1.
          Endif
          
          in  = (jp1-1)*ny + i
          is  = (jm1-1)*ny + i
          ie  = (j  -1)*ny + ip1
          iw  = (j  -1)*ny + im1
          
          ine = (jp1-1)*ny + ip1
          inw = (jp1-1)*ny + im1
          ise = (jm1-1)*ny + ip1
          isw = (jm1-1)*ny + im1

  !*** compute i-gradient

          If (.not. grid1_mask(ie)) then
             ie = n
             delew = 1.
          Endif
   
          If (.not. grid1_mask(iw)) then
             iw = n
             delew = 1.
          Endif

          gradi(n) = delew*(var1D(ie) - var1D(iw))

  !*** compute j-gradient
          
          If (.not. grid1_mask(in)) then
             in = n
             delns = 1.
          Endif
             
          If (.not. grid1_mask(is)) then
             is = n
             delns = 1.
          Endif
          
          gradj(n) = delns*(var1D(in) - var1D(is))                    

  !*** compute ij-gradient

          delew = 0.5
          
          If (jp1 == j .or. jm1 == j) then
             delns = 1.
          Else 
             delns = 0.5
          Endif
          
          If (.not. grid1_mask(ine)) then
             if (in /= n) then
                ine = in
                delew = 1.
             else if (ie /= n) then
                ine = ie
                inw = iw
                if (inw == n) delew = 1.
                delns = 1.
             else
                ine = n
                inw = iw
                delew = 1
                delns = 1
             endif
          Endif
          
          If (.not. grid1_mask(inw)) then
             if (in /= n) then
                inw = in
                delew = 1.
             else if (iw /= n) then
                inw = iw
                ine = ie
                if (ie == n) delew = 1.
                delns = 1.
             else
                inw = n
                ine = ie
                delew = 1.
                delns = 1.
             endif
          Endif
          
          deriv1(n) = delew*(var1D(ine)-var1D(inw))                    
          
          If (.not. grid1_mask(ise)) then
             if (is /= n) then
                ise = is
                delew = 1.
             else if (ie /= n) then
                ise = ie
                isw = iw
                if (isw == n) delew = 1.
                delns = 1.
             else
                ise = n
                isw = iw
                delew = 1.
                delns = 1.
             endif
          Endif
          
          If (.not. grid1_mask(isw)) then
             if (is /= n) then
                isw = is
                delew = 1.
             else if (iw /= n) then
                isw = iw
                ise = ie
                if (ie == n) delew = 1.
                delns = 1.
             else
                isw = n
                ise = ie
                delew = 1.
                delns = 1.
             endif
          Endif
          
          deriv2(n) = delew*(var1D(ise) - var1D(isw))
          gradij(n) = delns*(deriv1(n) - deriv2(n))
       ENDIF
    End do
    
    deallocate(deriv1,deriv2,grid1_mask)   
    allocate(var_interp1D(nxfin*nyfin))
    
    var_interp1D = 0.0
    num_links = size(dst_add)
    
    Do i = 1,num_links
       
       var_interp1D(dst_add(i)) = var_interp1D(dst_add(i))        +      &
            matrix(1,i)*var1D(src_add(i))   +      &
            matrix(2,i)*gradi(src_add(i))   +      &
            matrix(3,i)*gradj(src_add(i))   +      & 
            matrix(4,i)*gradij(src_add(i))
    End do
    
    deallocate(gradi,gradj,gradij,var1D)
    
    Call tab1Dto2D(var_interp1D,tabout,nyfin,nxfin)
    
    deallocate(var_interp1D)       
    
  end subroutine make_bicubic_remap
          
  
  !**************************************************************
  !   end subroutine make_bicubic_remap
  !**************************************************************          
  
end module outil
