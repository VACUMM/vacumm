% -- Programme validation MARS

% -- Module 1: colocalisation

% (suit le programme main_colocalisation.f90)
%
% Read namelist


% read_io.f90 -> r_obs.f90 et r_mars.f90



% 
% 
% % -----------------------------------------
% % 1) Read modelled fields and observations
% % read_io.f90 -> r_obs.f90 et r_mars.f90
% % -----------------------------------------
% Call read_io
% 
% % --------------------------------------------
% % 2) Interpolation obs/model - time dimension
% % --------------------------------------------
% Call time_interp
% 
% % --------------------------------------------
% % 2) Interpolation obs/model - space dimendion
% %
% % example: for SST -> lon/lat -> grid to grid
% % --------------------------------------------
% Call space_interp
% 
% 
% 
% 
