Release notes
#############

Version 3.3.0
=============

- Added cmap_lum, cma_sat and cmap_pastel keyword to plots
- Added saturate, desaturate, change_luminosity, change_saturation, pastelise
- Added cdtime validator to ConfigManager
- Added isempty to StatAccum
- Added interp_clim
- Added cylic mode support for extend1d
- Added support of intervals to get_xy
- Added espg support to cached_map
- Added x/ymargin support to minimap
- Added squarebox
- Added SimpleCloudKriger
- Added cellwidth support to bounds1d and meshcells
- Added dict support to scalebox
- Added cmocean colormaps support
- Added add_lightshading to plot
- Added fortran linear4dto1d
- Fixed xmlconfig attribute type checking
- Added zenodo DOI
- Renamed import of time_selector to filter_time_selector
- Fixed ncget_grid
- Fixed axis for 1d regridding of 1d arrays
- Fixed 360 wrap test in grid2xy
- Fixed tsel2slice
- Fixed template cloning in StatAccum
- Fixed dump/load of StatAccum
- Fixed resol_mask

Version 3.2.0
=============

- Added support for auto placement of text in add_place
- Added support for list of files to savefig
- Added suppor for dict to initialise childnodes in XmlConfig
- Added time split support to transect
- Added splitidx to get 1d splitting specs
- Added the add_grid method to Plot2D
- Added the x/ycorners options to add_grid function
- Added index support to ghhs_autores
- Added cmap+color validators to misc.config
- Added extraopts to ConfigManager
- Added get_quiverkey_value to plot
- Added x/y2db argument to Plot2D
- Added start_redirections and stop_redirections to log.Logger
- Fixed method testing in grid2xy and transect
- Fixed color in add_point
- Fixed masking in StepsNorm
- Fixed minute case in basic_auto_scaling
- Fixed generic2d masking and gaussian2d
- Fixed config read in grid module
- Fixed issue #2: verbose and notice fail for Logger subclasses

Version 3.1.1
=============

- Added redirection support to log.Logger.
- Added dstpts2line to interp.
- Added gen_binhelps extension to sphinx.
- Fixed format and date_format use in log.Logger.
- Fixed module members list in units.
- Fixed dstwgt2dto1dc_reduc.
- Fixed interp.mix2d for datarmor.
- Fixed StatAccum hist templates.
- Fixed gen_cmaps.
- Fixed plot_cmap and plot_cmaps.
- Fixed issue with Makefile.
- Fixed setup with CHANGES.

Version 3.1.0
=============

- Added support for mtype=None to variogram_fit.
- Added support for res="None" to create_map.
- Added errfunc support to kriging.
- Added cfgfilter to cfgmanager.
- Added proj param to basemap.get_proj.
- Added closing after showing in core_plot.
- Added autoscaling mode to ScalarMappable.get_levels with normal and degrees.
- Added merge_masks to merge masks of several variables.
- Added u, v, ubc, vb, speed, cdir, sigma*, *dens and renamed vol to cvol in cf.
- Added support for redirecting warnings, stdout and sterr to io.Logger
- Added mode support to dz2depth with edge, edge+ and middle.
- Added checkdir to make sure dir exists.
- Added julday converter.
- Added support for haversine distance to get_distances + krig integration.
- Added cyclic support to rainbow.
- Changed grid2xy to use get_distances.
- Renamed dmax to distmax in kriging.
- Improved support of julian days in atime.
- Improved date locators and formatters.
- Fixed cmap_br*.
- Fixed ignorecase in ncmatch_obj.
- Fixed some proj problems in misc.grid.
- Fixed inversions in kriging.
- Fixed validation of list in config.
- Fixed 360 deg problem for grid2xy.
- Fixed transect with 4D data
- Fixed scalar handling and masking in grid2xy
- Fixed format_var with variables with no axes specs.
- Fixed roundto in IterDates.
- Fixed no_norm issue.
- Fixed some standard names and grid locs in cf.

Version 3.0.0
=============

- Added font weight change for degrees in labels.
- Added standard_names to names for searching in cf.
- Added showvar.py to quickly display a netcdf variable.
- Added support for min+max+hist and restart to StatAccum.
- Added support for exact and block kriging to OCK.
- Added sill and range to linear variogram model in kriging.
- Added constraints to variogram model fit.
- Added color.discretize_cmap.
- Added Plot.add_water_mark.
- Added units.basic_proj.
- Added systematic cleaning to cache_map().
- Added [vacumm.misc.grid.basemap]max_cache_size config option.
- Added cellerr method to regrid1d.
- Added time arguments support if applicable to Plot.add_point().
- Added dstwgt method for fortran interpolators from gridded to random points.
- Added tuple support for time creation routines of atime.
- New regrid2d with tool and method keywords.
- Fixed range in hlitvs.
- Fixed mixed_layer_depth with kz.
- Fixed: default params in get_proj.
- Fixed names of module attributes which are now upper case.
- Fixed: list_forecast_files, Plot.add_lon/lat, _interp_.linept, Plot2D.fill.
- Fixed: ConfigManager.opt_parse.
- Removed sphinxfortran extension which is now a standalone vacumm project.

Version 2.5.4
=============

- Added "make safedoc" target.
- Fixed: english translations++.
- Fixed: missing test_plot_add_logo.py.
- Fixed: multifit+multiproc in kriging.
- Fixed: ConfigManager.arg_parse helps.
- Fixed: station_info import of oldnumeric.

Version 2.5.1
=============

- Changed: module level config files renamed to vacumm.cfg.
- Fixed: access to vacumm_nice_gfdl and vacumm_ssec colormaps.
- Fixed: Logger and specs for Profile.
- Fixed: add_logo.
- Fixed: removed dependency to pytz, which must now be installed
  to add time zone support to vacumm.

Version 2.5.0
==============

- Added: camp_nice_gfdl colormap.
- Added: Plot.add_annotation.
- Added: misc.plot.advanced.add_things tutorial.
- Fixed: gen_gallery.

Version 2.4.2
==============

- Added: misc.isempty.
- Fixed: cfg2rst, ConfigManager, StepsNorm.

Version 2.4.1
==============

- Upgraded: Logger.
- Added: docversions sphinx extension.
- Fixed: ConfigManager.opt_parse/arg_parse, Shapes, get_proj, get_xy,
  seawater import, are_good_units, Shapes.__init__/plot.

Version 2.4.0
==============

- Added: Added fp + th1p + some wind variables to cf.
- Added: add_arrow method to Plot2D.
- Added: add_map_places plot function.
- Improved: In curve2, an array can be passed to fill_between keyword.
- Fixed: ConfigManager, polygon_select, polygon_mask, coord2slice, sigma,
  tide.filters, StepsNorm, list_forecast_files, NEMO.

Version 2.3.1
=============

- Fortran regrid1d routines work directly with missing values.
- Unit tests save outputs in scripts directory.
- Fixed installation issue with setup.*.
- Fixed bugs: list_forecast_files, filter_selector, NEMO, coord2slice.

Version 2.3.0
=============

- Added the new CurvedInterpolator based on some fortran code
  primarily used for computing transects.
- New regrid1dnew that can regrid from a variable 1D axis to another
  variable 1D axis, like for instance from sigma to sigma coordinates.
  It will later replace regrid1d. Extrapolation in regrid1dnew is
  now available for all methods.
- Improvements for staggered grids in Dataset.
- minimap can now display background data instead of ocean color.
- cf: added wspd and wdir for wind.
- Smaller data samples.
- Better management of staggering in Dataset and arakawa (still experimental).
- Removed setup.cfg and added two templates, with a simple one and
  another one for OpenMP parallelisation.
- Fixed issues: vacumm config, sigma2depth, grid2xy, format_var,
  fortran_domain, etc.




