
module bicubic_interp
  
  use outil
  use declaration
  
  !************************************************************************
  ! 									*
  ! MODULE  BICUBIC INTERP						*
  ! 									*
  ! bicubic interpolation routines from SCRIP package			*		
  !									*
  !http://climate.lanl.gov/Software/SCRIP/				*
  !									*
  !Bicubic remapping							*
  !									*
  !************************************************************************
                 
  implicit none

  Integer :: grid1_size,grid2_size,grid1_rank, grid2_rank      
  Integer, dimension(:), pointer :: grid1_dims, grid2_dims  
  Logical, dimension(:), pointer :: grid1_mask,grid2_mask        
  Real(kind=rlg),dimension(:),pointer :: &
       grid1_center_lat,  &
       grid1_center_lon,  & 
       grid2_center_lat,  &
       grid2_center_lon,  &
       grid1_frac,        & 
       grid2_frac           
  Real(kind=rlg),dimension(:,:), pointer :: grid1_bound_box,grid2_bound_box   
  Integer, parameter :: num_srch_bins = 90  
  Integer,dimension(:,:),pointer :: bin_addr1,bin_addr2 
  Real(kind=rlg), dimension(:,:),pointer :: bin_lats,bin_lons       
  Real(kind=rlg), parameter :: zero   = 0.0,  &
       one    = 1.0,  &
       two    = 2.0,  &
       three  = 3.0,  &
       four   = 4.0,  &
       five   = 5.0,  & 
       half   = 0.5,  &
       quart  = 0.25, &
       bignum = 1.e+20, &
       tiny   = 1.e-14, &
       pi2    = two*pi, &
       pih    = half*pi        
      
  Real(kind=rlg), parameter :: deg2rad = pi/180.   
  Integer , parameter :: max_iter = 100   
  Real(kind=rlg), parameter :: converge = 1.e-10 
  Integer, parameter :: norm_opt_none    = 1 &
       ,norm_opt_dstarea = 2 &
       ,norm_opt_frcarea = 3
  
  Integer, parameter :: map_type_conserv  = 1 &
       ,map_type_bilinear = 2 &
       ,map_type_bicubic  = 3 &
       ,map_type_distwgt  = 4

  Integer :: max_links_map1  &  ! current size of link arrays
       ,num_links_map1  &  ! actual number of links for remapping
       ,max_links_map2  &  ! current size of link arrays
       ,num_links_map2  &  ! actual number of links for remapping
       ,num_maps        &  ! num of remappings for this grid pair
       ,num_wts         &  ! num of weights used in remapping
       ,map_type        &  ! identifier for remapping method
       ,norm_opt        &  ! option for normalization (conserv only)
       ,resize_increment ! default amount to increase array size
  
  Integer , dimension(:), pointer :: &
       bicubic_grid1_add_map1, &  ! grid1 address for each link in mapping 1
       bicubic_grid2_add_map1, &  ! grid2 address for each link in mapping 1
       grid1_add_map2, &  ! grid1 address for each link in mapping 2
       grid2_add_map2   ! grid2 address for each link in mapping 2
  
  Real(kind=rlg), dimension(:,:), pointer ::   &
       bicubic_wts_map1, &   ! map weights for each link (num_wts,max_links)
       wts_map2 => NULL()    ! map weights for each link (num_wts,max_links)

  !***********************************************************************

      contains

  !***********************************************************************
        subroutine get_remap_bicub(grid1_lat,grid2_lat,grid1_lon,grid2_lon,mask,maskdst, &
             remap_matrix,source_add,destination_add)
          
  !----------------------------------------------------------------------- 
  !
  !     this routine computes the weights for a bicubic interpolation.
  ! 
  !-----------------------------------------------------------------------
  
          Real(kind=rlg),dimension(:,:),pointer :: grid1_lat,grid2_lat,grid1_lon,grid2_lon
          Logical,dimension(:,:) :: mask,maskdst
          
          Integer,dimension(:),pointer :: source_add,destination_add
          Real(kind=rlg),dimension(:,:),pointer :: remap_matrix 

  !-----------------------------------------------------------------------
  !
  !     local variables
  !
  !-----------------------------------------------------------------------
          
          Integer :: n,icount,dst_add,iter,nx,ny,i,j,ip1,e_add,jp1,n_add,ne_add,nele
          Real(kind=rlg), dimension(4) :: tmp_lats, tmp_lons
          Integer, dimension(4) :: src_add
          
          Real(kind=rlg), dimension(4)  :: src_lats,src_lons
          Integer lastsrc_add
          Real(kind=rlg), dimension(4,4)  :: wgts
          Real(kind=rlg) :: dlat
          Real(kind=rlg) :: plat, plon,iguess, jguess,deli, delj, &
               dth1, dth2, dth3,dph1, dph2, dph3,dthp, dphp,          &
               mat1, mat2, mat3, mat4,determinant,sum_wgts                     
          
  !-----------------------------------------------------------------------
  !
  !     compute mappings from grid1 to grid2
  !
  !-----------------------------------------------------------------------

          allocate(grid1_dims(2),grid2_dims(2))
          grid1_dims(1) = size(grid1_lat,2)
          grid1_dims(2) = size(grid1_lat,1)
          grid2_dims(1) = size(grid2_lat,2)
          grid2_dims(2) = size(grid2_lat,1)
          grid1_size = size(grid1_lat,2) * size(grid1_lat,1)
          grid2_size = size(grid2_lat,2) * size(grid2_lat,1) 
      
          allocate( grid2_mask(grid2_size),         &
               grid1_bound_box (4,grid1_size), &
               grid2_bound_box (4,grid2_size), &
               grid1_frac      (grid1_size),   &
               grid2_frac      (grid2_size))
          
          grid1_frac = zero
          grid2_frac = zero
          
  !
  ! 2D array -> 1D array
  !
      allocate(grid1_center_lat(size(grid1_lat,1)*size(grid1_lat,2)))
      call tab2Dto1D(grid1_lat,grid1_center_lat)
      
      allocate(grid1_center_lon(size(grid1_lon,1)*size(grid1_lon,2)))
      call tab2Dto1D(grid1_lon,grid1_center_lon)
      
      allocate(grid2_center_lat(size(grid2_lat,1)*size(grid2_lat,2)))
      call tab2Dto1D(grid2_lat,grid2_center_lat)
      
      allocate(grid2_center_lon(size(grid2_lon,1)*size(grid2_lon,2)))      
      call tab2Dto1D(grid2_lon,grid2_center_lon) 
      
      allocate(grid1_mask(size(grid1_lat,1)*size(grid1_lat,2)))
      call logtab2Dto1D(mask,grid1_mask)
      
      call logtab2Dto1D(maskdst,grid2_mask)

      
  !      Write(*,*) ,'grid1_mask = ',grid1_mask 
              

  ! degrees to radian

      grid1_center_lat = grid1_center_lat*deg2rad
      grid1_center_lon = grid1_center_lon*deg2rad
      grid2_center_lat = grid2_center_lat*deg2rad
      grid2_center_lon = grid2_center_lon*deg2rad
      
  !-----------------------------------------------------------------------
  !     convert longitudes to 0,2pi interval
  !-----------------------------------------------------------------------
      
      where (grid1_center_lon .gt. pi2)  grid1_center_lon =       &
           grid1_center_lon - pi2
      where (grid1_center_lon .lt. zero) grid1_center_lon =       &
           grid1_center_lon + pi2
      where (grid2_center_lon .gt. pi2)  grid2_center_lon =       &
           grid2_center_lon - pi2
      where (grid2_center_lon .lt. zero) grid2_center_lon =       &
           grid2_center_lon + pi2

  !-----------------------------------------------------------------------
  !
  !     make sure input latitude range is within the machine values
  !     for +/- pi/2 
  !
  !-----------------------------------------------------------------------

      where (grid1_center_lat >  pih) grid1_center_lat =  pih
      where (grid1_center_lat < -pih) grid1_center_lat = -pih
      where (grid2_center_lat >  pih) grid2_center_lat =  pih
      where (grid2_center_lat < -pih) grid2_center_lat = -pih
      
  !-----------------------------------------------------------------------
  !
  !     compute bounding boxes for restricting future grid searches
  !
  !-----------------------------------------------------------------------

      nx = grid1_dims(1)
      ny = grid1_dims(2)
      
      do n=1,grid1_size
         
  !*** find N,S and NE points to this grid point
         
         j = (n - 1)/nx +1
         i = n - (j-1)*nx
         
         if (i < nx) then
            ip1 = i + 1
         else
  !*** assume cyclic
            ip1 = 1
  !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid1_center_lat(e_add) -     &
                 grid1_center_lat(n   )) > pih) then
               ip1 = i
            endif
            ip1=nx
         endif
         
         if (j < ny) then
            jp1 = j+1
         else
  !*** assume cyclic
            jp1 = 1
  !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid1_center_lat(n_add) -             &
                 grid1_center_lat(n   )) > pih) then
               jp1 = j
            endif
            jp1=ny
         endif
         
         n_add = (jp1 - 1)*nx + i
         e_add = (j - 1)*nx + ip1
         ne_add = (jp1 - 1)*nx + ip1
         
  !*** find N,S and NE lat/lon coords and check bounding box
         
         tmp_lats(1) = grid1_center_lat(n)
         tmp_lats(2) = grid1_center_lat(e_add)
         tmp_lats(3) = grid1_center_lat(ne_add)
         tmp_lats(4) = grid1_center_lat(n_add)
         
         tmp_lons(1) = grid1_center_lon(n)
         tmp_lons(2) = grid1_center_lon(e_add)
         tmp_lons(3) = grid1_center_lon(ne_add)
         tmp_lons(4) = grid1_center_lon(n_add)
         
         grid1_bound_box(1,n) = minval(tmp_lats)
         grid1_bound_box(2,n) = maxval(tmp_lats)
         
         grid1_bound_box(3,n) = minval(tmp_lons)
         grid1_bound_box(4,n) = maxval(tmp_lons)
        
      end do
      
      nx = grid2_dims(1)
      ny = grid2_dims(2)
      
      do n=1,grid2_size
         
  !*** find N,S and NE points to this grid point

         j = (n - 1)/nx +1
         i = n - (j-1)*nx
         
         if (i < nx) then
            ip1 = i + 1
         else
  !*** assume cyclic
            ip1 = 1
  !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid2_center_lat(e_add) -  &
                 grid2_center_lat(n   )) > pih) then
               ip1 = i
            endif
         endif
         
         if (j < ny) then
            jp1 = j+1
         else
  !*** assume cyclic
            jp1 = 1
  !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid2_center_lat(n_add) -  &
                 grid2_center_lat(n   )) > pih) then
               jp1 = j
            endif
         endif
         
         n_add = (jp1 - 1)*nx + i
         e_add = (j - 1)*nx + ip1
         ne_add = (jp1 - 1)*nx + ip1
         
  !*** find N,S and NE lat/lon coords and check bounding box

         tmp_lats(1) = grid2_center_lat(n)
         tmp_lats(2) = grid2_center_lat(e_add)
         tmp_lats(3) = grid2_center_lat(ne_add)
         tmp_lats(4) = grid2_center_lat(n_add)
         
         tmp_lons(1) = grid2_center_lon(n)
         tmp_lons(2) = grid2_center_lon(e_add)
         tmp_lons(3) = grid2_center_lon(ne_add)
         tmp_lons(4) = grid2_center_lon(n_add)
         
         grid2_bound_box(1,n) = minval(tmp_lats)
         grid2_bound_box(2,n) = maxval(tmp_lats)
         grid2_bound_box(3,n) = minval(tmp_lons)
         grid2_bound_box(4,n) = maxval(tmp_lons)
      end do

      where (abs(grid1_bound_box(4,:) - grid1_bound_box(3,:)) > pi)
         grid1_bound_box(3,:) = zero
         grid1_bound_box(4,:) = pi2
      end where
      
      where (abs(grid2_bound_box(4,:) - grid2_bound_box(3,:)) > pi)
         grid2_bound_box(3,:) = zero
         grid2_bound_box(4,:) = pi2
      end where
      
  !***
  !*** try to check for cells that overlap poles
  !***
      
      where (grid1_center_lat > grid1_bound_box(2,:)) &
           grid1_bound_box(2,:) = pih
      
      where (grid1_center_lat < grid1_bound_box(1,:)) &
           grid1_bound_box(1,:) = -pih
      
      where (grid2_center_lat > grid2_bound_box(2,:)) &
           grid2_bound_box(2,:) = pih
      
      where (grid2_center_lat < grid2_bound_box(1,:)) &
           grid2_bound_box(1,:) = -pih
      
  !-----------------------------------------------------------------------
  !     set up and assign address ranges to search bins in order to 
  !     further restrict later searches
  !-----------------------------------------------------------------------
      
      allocate(bin_addr1(2,num_srch_bins))
      allocate(bin_addr2(2,num_srch_bins))
      allocate(bin_lats (2,num_srch_bins))
      allocate(bin_lons (2,num_srch_bins))
      
      dlat = pi/num_srch_bins
      
      do n=1,num_srch_bins
         bin_lats(1,n) = (n-1)*dlat - pih
         bin_lats(2,n) =     n*dlat - pih
         bin_lons(1,n) = zero
         bin_lons(2,n) = pi2
         bin_addr1(1,n) = grid1_size + 1
         bin_addr1(2,n) = 0
         bin_addr2(1,n) = grid2_size + 1
         bin_addr2(2,n) = 0
      end do
      
      do nele=1,grid1_size
         do n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and.   &
                 grid1_bound_box(2,nele) >= bin_lats(1,n)) then
               bin_addr1(1,n) = min(nele,bin_addr1(1,n))
               bin_addr1(2,n) = max(nele,bin_addr1(2,n))
            endif
         end do
      end do
      
      do nele=1,grid2_size
         do n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and.    &
                 grid2_bound_box(2,nele) >= bin_lats(1,n)) then
               bin_addr2(1,n) = min(nele,bin_addr2(1,n))
               bin_addr2(2,n) = max(nele,bin_addr2(2,n))
            endif
         end do
      end do
      
      Call init_remap_vars    
      
  !***
  !*** loop over destination grid 
  !***
      
      !grid2_mask = .true.

      lastsrc_add=1
    
      Write(*,*) 'Bicubic remapping weights computation'
 
      grid_loop1: DO dst_add = 1, grid2_size
         
         if (.not. grid2_mask(dst_add)) then
            cycle grid_loop1
         end if
         
         
         plat = grid2_center_lat(dst_add)
         plon = grid2_center_lon(dst_add)         
          
                    
         Call grid_search_bicub(src_add, src_lats, src_lons,     &
              plat, plon, grid1_dims,              &
              grid1_center_lat, grid1_center_lon,  &
              grid1_bound_box, bin_addr1, lastsrc_add)
              !grid1_bound_box, bin_addr1, bin_addr2,lastsrc_add)
            
         If (src_add(1) > 0) then
            Do n=1,4
               if (.not. grid1_mask(src_add(n))) src_add(1) = 0
            End Do
         Endif
         
  !-----------------------------------------------------------------------
  !
  !       if point found, find local i,j coordinates for weights
  !
  !-----------------------------------------------------------------------
         
         IF (src_add(1) > 0) then
            
            grid2_frac(dst_add) = one
            
            dth1 = src_lats(2) - src_lats(1)
            dth2 = src_lats(4) - src_lats(1)
            dth3 = src_lats(3) - src_lats(2) - dth2
            
            dph1 = src_lons(2) - src_lons(1)
            dph2 = src_lons(4) - src_lons(1)
            dph3 = src_lons(3) - src_lons(2)
            
            if (dph1 >  three*pih) dph1 = dph1 - pi2
            if (dph2 >  three*pih) dph2 = dph2 - pi2
            if (dph3 >  three*pih) dph3 = dph3 - pi2
            if (dph1 < -three*pih) dph1 = dph1 + pi2
            if (dph2 < -three*pih) dph2 = dph2 + pi2
            if (dph3 < -three*pih) dph3 = dph3 + pi2
            
            dph3 = dph3 - dph2
            iguess = half
            jguess = half
            
            Do iter=1,max_iter
               
               dthp = plat-src_lats(1)-dth1*iguess- &
                    dth2*jguess-dth3*iguess*jguess
               dphp = plon - src_lons(1)
               
               if (dphp >  three*pih) dphp = dphp - pi2
               if (dphp < -three*pih) dphp = dphp + pi2
               
               dphp = dphp - dph1*iguess - dph2*jguess - & 
                    dph3*iguess*jguess
               
               mat1 = dth1 + dth3*jguess
               mat2 = dth2 + dth3*iguess
               mat3 = dph1 + dph3*jguess
               mat4 = dph2 + dph3*iguess
               
               determinant = mat1*mat4 - mat2*mat3
               
               deli = (dthp*mat4 - mat2*dphp)/determinant
               delj = (mat1*dphp - dthp*mat3)/determinant
               
               if (abs(deli) < converge .and.            & 
                    abs(delj) < converge) exit
               
               iguess = iguess + deli
               jguess = jguess + delj
               
            End Do
            
            If (iter <= max_iter) then
               
               wgts(1,1) = (one - jguess**2*(three-two*jguess))*   &
                    (one - iguess**2*(three-two*iguess))
               
               wgts(1,2) = (one - jguess**2*(three-two*jguess))*   &
                    iguess**2*(three-two*iguess)
               
               wgts(1,3) = jguess**2*(three-two*jguess)*           &
                    iguess**2*(three-two*iguess)
               
               wgts(1,4) = jguess**2*(three-two*jguess)*           &
                    (one - iguess**2*(three-two*iguess))
               
               wgts(2,1) = (one - jguess**2*(three-two*jguess))*   &
                    iguess*(iguess-one)**2
               
               wgts(2,2) = (one - jguess**2*(three-two*jguess))*   &
                                        iguess**2*(iguess-one)
               
               wgts(2,3) = jguess**2*(three-two*jguess)*           &
                    iguess**2*(iguess-one)
               
               wgts(2,4) =  jguess**2*(three-two*jguess)*          &
                    iguess*(iguess-one)**2
               
               wgts(3,1) =  jguess*(jguess-one)**2*                &
                    (one - iguess**2*(three-two*iguess))
               
               wgts(3,2) =  jguess*(jguess-one)**2*                &
                    iguess**2*(three-two*iguess)                            
               
               wgts(3,3) =  jguess**2*(jguess-one)*                &
                    iguess**2*(three-two*iguess)
               
               wgts(3,4) =  jguess**2*(jguess-one)*                &
                    (one - iguess**2*(three-two*iguess))
               
               wgts(4,1) =  iguess*(iguess-one)**2*                &
                    jguess*(jguess-one)**2
               
               wgts(4,2) =  iguess**2*(iguess-one)*                &
                    jguess*(jguess-one)**2
               
               wgts(4,3) =  iguess**2*(iguess-one)*                &
                    jguess**2*(jguess-one)
               
               wgts(4,4) =  iguess*(iguess-one)**2*                &
                    jguess**2*(jguess-one)
               
               call store_link_bicub(dst_add, src_add, wgts)
               
            Else
               stop 'Iteration for i,j exceed max iteration count'
            Endif
            
         ELSE IF (src_add(1) < 0) then
            
            src_add = abs(src_add)
            
            icount = 0
            
            Do n=1,4
               if (grid1_mask(src_add(n))) then
                  icount = icount + 1
               else
                  src_lats(n) = zero
               endif
            End Do
            
            If (icount > 0) then
               
               sum_wgts = sum(src_lats)
               wgts(1,1) = src_lats(1)/sum_wgts
               wgts(1,2) = src_lats(2)/sum_wgts
               wgts(1,3) = src_lats(3)/sum_wgts
               wgts(1,4) = src_lats(4)/sum_wgts
               wgts(2:4,:) = zero
     !          wgts(2,1) = src_lats(1)/sum_wgts
     !          wgts(2,2) = src_lats(2)/sum_wgts
     !          wgts(2,3) = src_lats(3)/sum_wgts
     !          wgts(2,4) = src_lats(4)/sum_wgts
     !          
     !          wgts(3,1) = src_lats(1)/sum_wgts
     !          wgts(3,2) = src_lats(2)/sum_wgts
     !          wgts(3,3) = src_lats(3)/sum_wgts
     !          wgts(3,4) = src_lats(4)/sum_wgts
     !          wgts(4,1) = src_lats(1)/sum_wgts
     !          wgts(4,2) = src_lats(2)/sum_wgts
     !          wgts(4,3) = src_lats(3)/sum_wgts
     !          wgts(4,4) = src_lats(4)/sum_wgts
               Call store_link_bicub(dst_add, src_add, wgts)
               
            Endif
            
         ENDIF
         
      END DO grid_loop1
      
      Allocate(remap_matrix(size(bicubic_wts_map1,1),size(bicubic_wts_map1,2)), &
           source_add(num_links_map1),       &
           destination_add(num_links_map1))
      
      do j = 1,size(bicubic_wts_map1,2)
         do i = 1,size(bicubic_wts_map1,1)
            
            remap_matrix(i,j) = bicubic_wts_map1(i,j)
            
         end do
      end do

  !print*,remap_matrix(1,1:30)
      

      source_add(:) = bicubic_grid1_add_map1(1:num_links_map1)
      destination_add(:) = bicubic_grid2_add_map1(1:num_links_map1)
           
      Where(destination_add == 0)
         destination_add = 1
      End Where
      
      Where(source_add == 0)
         source_add = 1
      End Where
         
      deallocate(grid1_bound_box,grid2_bound_box,grid1_center_lat,grid1_center_lon)
      deallocate(grid2_center_lat,grid2_center_lon,bicubic_grid1_add_map1,bicubic_wts_map1)
      deallocate(grid1_frac,grid2_frac,grid1_dims,grid2_dims,grid2_mask)
      deallocate(bin_addr1,bin_addr2,bin_lats,bin_lons)
      deallocate(grid1_mask,bicubic_grid2_add_map1)	

  !-----------------------------------------------------------------------
      
    end subroutine get_remap_bicub  

  !***********************************************************************
    
  !***********************************************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************************************************************************
  !   SUBROUTINE GRID_SEARCH_BILIN
  !************************************************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine grid_search_bicub(src_add, src_lats, src_lons,   &
         plat, plon, src_grid_dims,      &
         src_center_lat, src_center_lon, & 
         src_grid_bound_box,             &
         src_bin_add, lastsrc_add)
         !src_bin_add, dst_bin_add,lastsrc_add)

  !----------------------------------------------------------------------- 
  !
  !     this routine finds the location of the search point plat, plon
  !     in the source grid and returns the corners needed for a bilinear
  !     interpolation.
  !
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !     output variables
  !-----------------------------------------------------------------------
  !
  ! address of each corner point enclosing P
  !
      integer,dimension(4) :: src_add  
      real(kind=rlg),dimension(4) :: src_lats,src_lons  
  !      
  !-----------------------------------------------------------------------
  !     input variables
  !-----------------------------------------------------------------------
 
  ! latitude, longitude of the search point

      real(kind=rlg), intent(in) :: plat,plon   

  ! size of each src grid dimension

      integer, dimension(2), intent(in) :: src_grid_dims  

  ! latitude, longitude of each src grid center

      real(kind=rlg), dimension(:), intent(in) :: src_center_lat,src_center_lon  

  ! bound box for source grid

      real(kind=rlg), dimension(:,:), intent(in) :: src_grid_bound_box 
 
  ! latitude bins for restricting searches

      integer, dimension(:,:), intent(in) ::src_bin_add 
      !integer, dimension(:,:), intent(in) ::dst_bin_add 
      
      integer,optional :: lastsrc_add
      integer :: loopsrc,l1,l2
 
  !-----------------------------------------------------------------------
  !     local variables
  !-----------------------------------------------------------------------

      integer :: n,next_n,srch_add,nx, ny,min_add, max_add,  &
           i, j, jp1, ip1, n_add, e_add, ne_add
      
      
      real(kind=rlg) ::  vec1_lat, vec1_lon,vec2_lat, vec2_lon, cross_product,  &
           cross_product_last,coslat_dst, sinlat_dst, coslon_dst, &
           sinlon_dst,dist_min, distance 
      
  !-----------------------------------------------------------------------
  !     restrict search first using bins
  !-----------------------------------------------------------------------
      src_add = 0
      min_add = size(src_center_lat)
      max_add = 1
      do n=1,num_srch_bins
         if (plat >= bin_lats(1,n) .and. plat <= bin_lats(2,n) .and. &
              plon >= bin_lons(1,n) .and. plon <= bin_lons(2,n)) then
            min_add = min(min_add, src_bin_add(1,n))
            max_add = max(max_add, src_bin_add(2,n))
         endif
      end do
      
  !-----------------------------------------------------------------------
  !     now perform a more detailed search 
  !-----------------------------------------------------------------------
      
      nx = src_grid_dims(1)
      ny = src_grid_dims(2)
      
      loopsrc=0
      do while (loopsrc <= max_add)
         
         
         l1=max(min_add,lastsrc_add-loopsrc)
         l2=min(max_add,lastsrc_add+loopsrc)      
         
         loopsrc = loopsrc+1
         
         srch_loop: do srch_add = l1,l2,max(l2-l1,1)
            
  !*** first check bounding box
        
            if (plat <= src_grid_bound_box(2,srch_add) .and. & 
                 plat >= src_grid_bound_box(1,srch_add) .and.  &
                 plon <= src_grid_bound_box(4,srch_add) .and.  &
                 plon >= src_grid_bound_box(3,srch_add)) then
               
  !***
  !*** we are within bounding box so get really serious
  !***
  !*** determine neighbor addresses
               
               j = (srch_add - 1)/nx +1
               i = srch_add - (j-1)*nx
               
               if (i < nx) then
                  ip1 = i + 1
               else
                  ip1 = 1
               endif
               
               if (j < ny) then
                  jp1 = j+1
               else
                  jp1 = 1
               endif
               
               n_add = (jp1 - 1)*nx + i
               e_add = (j - 1)*nx + ip1
               ne_add = (jp1 - 1)*nx + ip1
               
               src_lats(1) = src_center_lat(srch_add)
               src_lats(2) = src_center_lat(e_add)
               src_lats(3) = src_center_lat(ne_add)
               src_lats(4) = src_center_lat(n_add)
               
               src_lons(1) = src_center_lon(srch_add)
               src_lons(2) = src_center_lon(e_add)
               src_lons(3) = src_center_lon(ne_add)
               src_lons(4) = src_center_lon(n_add)
  !
  !***
  !*** for consistency, we must make sure all lons are in
  !*** same 2pi interval
  !***
  !
               vec1_lon = src_lons(1) - plon
               if (vec1_lon >  pi) then
                  src_lons(1) = src_lons(1) - pi2
               else if (vec1_lon < -pi) then
                  src_lons(1) = src_lons(1) + pi2
               endif
               do n=2,4
                  vec1_lon = src_lons(n) - src_lons(1)
                  if (vec1_lon >  pi) then
                     src_lons(n) = src_lons(n) - pi2
                  else if (vec1_lon < -pi) then
                     src_lons(n) = src_lons(n) + pi2
                  endif
               end do

               corner_loop: do n=1,4
                  next_n = MOD(n,4) + 1
  !***
  !*** here we take the cross product of the vector making 
  !*** up each box side with the vector formed by the vertex
  !*** and search point.  if all the cross products are 
  !*** positive, the point is contained in the box.
  !***
                  vec1_lat = src_lats(next_n) - src_lats(n)
                  vec1_lon = src_lons(next_n) - src_lons(n)
                  vec2_lat = plat - src_lats(n)
                  vec2_lon = plon - src_lons(n)
  !***
  !*** check for 0,2pi crossings
  !***
                  if (vec1_lon >  three*pih) then
                     vec1_lon = vec1_lon - pi2
                  else if (vec1_lon < -three*pih) then
                     vec1_lon = vec1_lon + pi2
                  endif
                  if (vec2_lon >  three*pih) then
                     vec2_lon = vec2_lon - pi2
                  else if (vec2_lon < -three*pih) then
                     vec2_lon = vec2_lon + pi2
                  endif
                  
                  cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat
  !
  !***
  !*** if cross product is less than zero, this cell
  !*** doesn't work
  !***
                  if (n == 1) cross_product_last = cross_product
                  if (cross_product*cross_product_last < zero) &
                       exit corner_loop
                  cross_product_last = cross_product
                  
               end do corner_loop
  !***
  !*** if cross products all same sign, we found the location
  !***
               if (n > 4) then
                  src_add(1) = srch_add
                  src_add(2) = e_add
                  src_add(3) = ne_add
                  src_add(4) = n_add
                  
                  lastsrc_add = srch_add
                  return
               endif
  !***
  !*** otherwise move on to next cell
  !***
            endif !bounding box check
      end do srch_loop
         
   enddo

  !***
  !*** if no cell found, point is likely either in a box that
  !*** straddles either pole or is outside the grid.  fall back
  !*** to a distance-weighted average of the four closest
  !*** points.  go ahead and compute weights here, but store
  !*** in src_lats and return -add to prevent the parent
  !*** routine from computing bilinear weights
  !***
  !print *,'Could not find location for ',plat,plon
  !print *,'Using nearest-neighbor average for this point'

   coslat_dst = cos(plat)
   sinlat_dst = sin(plat)
   coslon_dst = cos(plon)
   sinlon_dst = sin(plon)

   dist_min = bignum
   src_lats = bignum
   
   do srch_add = min_add,max_add
      distance = acos(coslat_dst*cos(src_center_lat(srch_add))*   &
           (coslon_dst*cos(src_center_lon(srch_add)) +   &
           sinlon_dst*sin(src_center_lon(srch_add)))+   &
           sinlat_dst*sin(src_center_lat(srch_add)))
      
      if (distance < dist_min) then
         sort_loop: do n=1,4
            if (distance < src_lats(n)) then
               do i=4,n+1,-1
                  src_add (i) = src_add (i-1)
                  src_lats(i) = src_lats(i-1)
               end do
               src_add (n) = -srch_add
               src_lats(n) = distance
               dist_min = src_lats(4)
               exit sort_loop
            endif
         end do sort_loop
      endif
   end do
   
   src_lons = one/(src_lats + tiny)
   distance = sum(src_lons)
   src_lats = src_lons/distance

  !-----------------------------------------------------------------------

 end subroutine grid_search_bicub

  !-----------------------------------------------------------------------

 subroutine store_link_bicub(dst_add, src_add, weights)

  !----------------------------------------------------------------------- 
  !
  !     this routine stores the address and weight for four links 
  !     associated with one destination point in the appropriate address 
  !     and weight arrays and resizes those arrays if necessary.
  !
  !-----------------------------------------------------------------------
   
  !-----------------------------------------------------------------------
  !
  !     input variables
  !
  !-----------------------------------------------------------------------
   
   Integer :: dst_add
   Integer, dimension(4) :: src_add
   Real(kind=rlg), dimension(4,4) :: weights 
   
  !-----------------------------------------------------------------------
  !
  !     local variables
  !
  !-----------------------------------------------------------------------
   
   Integer :: n,num_links_old
   
  !-----------------------------------------------------------------------
  !
  !     increment number of links and check to see if remap arrays need
  !     to be increased to accomodate the new link.  then store the
  !     link.
  !
  !-----------------------------------------------------------------------
   
   
   num_links_old  = num_links_map1
   num_links_map1 = num_links_old + 4
   
   if (num_links_map1 > max_links_map1)    & 
        call resize_remap_vars(resize_increment)
   
   do n=1,4
      bicubic_grid1_add_map1(num_links_old+n) = src_add(n)
      bicubic_grid2_add_map1(num_links_old+n) = dst_add
      bicubic_wts_map1    (:,num_links_old+n) = weights(:,n)
   end do
   
  !----------------------------------------------------------------------- 
   
 end subroutine store_link_bicub
 
  !***********************************************************************
  !
  !
  !************************************************************************
  !   SUBROUTINE INIT_REMAP_VAR
  !************************************************************************

 subroutine init_remap_vars

   num_wts = 4
   
   num_links_map1 = 0
   max_links_map1 = 4*grid2_size
   if (num_maps > 1) then
      num_links_map2 = 0
      max_links_map1 = max(4*grid1_size,4*grid2_size)
      max_links_map2 = max_links_map1
   endif
   
   resize_increment = 0.1*max(grid1_size,grid2_size)
   
  !-----------------------------------------------------------------------
  !
  !     allocate address and weight arrays for mapping 1
  !   
  !-----------------------------------------------------------------------
  
   allocate (bicubic_grid1_add_map1(max_links_map1),  &
        bicubic_grid2_add_map1(max_links_map1),  &
        bicubic_wts_map1(num_wts, max_links_map1))
   
  !-----------------------------------------------------------------------
  !
  !     allocate address and weight arrays for mapping 2 if necessary 
  !   
  !-----------------------------------------------------------------------
   
   if (num_maps > 1) then
      allocate (grid1_add_map2(max_links_map2),    &
           grid2_add_map2(max_links_map2),    &
           wts_map2(num_wts, max_links_map2))
   endif
   
  !-----------------------------------------------------------------------
   
 end subroutine init_remap_vars

  !************************************************************************
  !************************************************************************
  !   SUBROUTINE RESIZE_REMAP_VAR
  !************************************************************************

 subroutine resize_remap_vars(increment)
  !-----------------------------------------------------------------------
  !     this routine resizes remapping arrays by increasing(decreasing)
  !     the max_links by increment
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !     input variables
  !-----------------------------------------------------------------------
   
   integer :: increment
   
  !-----------------------------------------------------------------------
  !      local variables
  !-----------------------------------------------------------------------
   
   integer ::    &
        mxlinks   ! size of link arrays
   
   integer, dimension(:), pointer ::    &
        add1_tmp,   & ! temp array for resizing address arrays
        add2_tmp  ! temp array for resizing address arrays
  !
  ! temp array for resizing weight arrays
  !
   real(kind=rlg), dimension(:,:), pointer :: wts_tmp   
  !
  !-----------------------------------------------------------------------
  !***
  !*** allocate temporaries to hold original values
  !***
  print*,'passage rezize' 
   mxlinks = size(bicubic_grid1_add_map1)
   allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), &
        wts_tmp(num_wts,mxlinks))
   
   add1_tmp = bicubic_grid1_add_map1
   add2_tmp = bicubic_grid2_add_map1
   wts_tmp  = bicubic_wts_map1
   
  !***
  !*** deallocate originals and increment max_links then
  !*** reallocate arrays at new size
  !***

   deallocate (bicubic_grid1_add_map1, bicubic_grid2_add_map1, bicubic_wts_map1)
   max_links_map1 = mxlinks + increment
   allocate (bicubic_grid1_add_map1(max_links_map1),    &
        bicubic_grid2_add_map1(max_links_map1),    &
        bicubic_wts_map1(num_wts,max_links_map1))
  !***
  !*** restore original values from temp arrays and
  !*** deallocate temps
  !***
   mxlinks = min(mxlinks, max_links_map1)
   bicubic_grid1_add_map1(1:mxlinks) = add1_tmp (1:mxlinks)
   bicubic_grid2_add_map1(1:mxlinks) = add2_tmp (1:mxlinks)
   bicubic_wts_map1    (:,1:mxlinks) = wts_tmp(:,1:mxlinks)
   deallocate(add1_tmp, add2_tmp, wts_tmp)

  !-----------------------------------------------------------------------
  !
 end subroutine resize_remap_vars
  !
  !************************************************************************
  !
 
end module bicubic_interp
