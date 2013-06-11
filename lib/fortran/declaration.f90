  !************************************************************
  ! Module containing variables declarations and main parameters
  ! used in all subroutines
  ! Created by M.Caillaud 01/02/09
  !*************************************************************

module declaration

  implicit none

  integer,parameter :: rsh=8,rlg=8
  integer,parameter :: lchain1=100
  integer :: imin=0,jmin=0
  integer :: rank
  real(kind=rlg), parameter ::pi=3.141592653589793238462643383279_rlg
  real(kind=rlg),parameter :: rad=180/pi
  real(kind=rsh) :: landvalue
  real(kind=rsh),parameter :: valmanq=-999.
  real(kind=rsh),parameter :: flag=500
  real(kind=rsh),parameter :: hmin=-15.

  ! ** NETCDF PARAMETERS **
  
  character(len=lchain1),parameter :: dimlon_name="longitude",    &
       dimlat_name="latitude",    &
       lon_name="longitude", &
       lonU_name="longitude_U", &
       lonV_name="longitude_V", &
       lat_name="latitude",  &
       latU_name="latitude_U",  &
       latV_name="latitude_V",  &
       H0_name="H0",         &
       HX_name="HX" ,        &
       HY_name="HY" ,         &
       meanlevel_name='nivmoy'
  integer :: status,ncid,h0_id,hx_id,hy_id,lon_id,lat_id,dimlon_id,dimlat_id,meanlevel_id
  integer :: lonU_id,latV_id
  
  ! ** namelist paths

  logical :: interp_sondes2grid,interp_grid2grid,lsmooth,lconnect,bmg,bmg_mask,l_closed_line
  character(len=lchain1) :: head_path,     &
       data_catalog,&
       data_path,   &
       nivmoypath    
  character(len=lchain1) :: output_file
  character(len=lchain1) :: coastfile
  character(len=4) ::  mask_method
  character(len=lchain1),parameter ::namelistname='./namelist'
  character(len=lchain1) :: connect_file1,connect_file2
  character(len=lchain1) ::  smooth_file
  logical :: south_border,north_border,east_border,west_border
  real(kind=rsh) :: rmax
  integer :: L,nbx,nby 
  logical :: l_bathy_threshold
  real(kind=rlg) :: bathy_threshold
  
end module declaration
