#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, inspect, os, sys
import cdms2, MV2
import vacumm.misc.bases as bases
import vacumm.misc.io as io
import vacumm.misc.grid.masking as masking

class Mask(bases.Object):

	inpfile = None
	outfile = None
	overwrite = None
	included_variables = None
	excluded_variables = None
	masked_only = None
	netcdf3 = None
	netcdf4level = 3
	
	argspec = inspect.getargspec(masking.polygon_mask)
	polygon_mask_defaults = dict(zip(argspec.args[-len(argspec.defaults):], argspec.defaults))
	del argspec
	
	resolution = 'i'
	resolutions = ('c','l','i','f','h')
	mode = polygon_mask_defaults.get('mode', 'i')
	modes = ('inside', 'intersect')
	thresholds = polygon_mask_defaults.get('thresholds', (0.5, 0.75))
	
	reverse = None
	
	@classmethod
	def get_argparser(cls):
		parser = argparse.ArgumentParser(description='Mask variables of a NetCDF file')
		parser.add_argument('-i', '--inpfile', required=True, help='NetCDF input file')
		parser.add_argument('-o', '--outfile', required=True, help='NetCDF output file')
		parser.add_argument('-O', '--overwrite', action='store_true', help='Overwrite output file if it already exists')
		vars_group = parser.add_mutually_exclusive_group()
		vars_group.add_argument('-I', '--include', dest='included_variables', action='append', help='Process only this variable (repeatable option)')
		vars_group.add_argument('-E', '--exclude', dest='excluded_variables', action='append', help='Process all except this variable (repeatable option)')
		parser.add_argument('-M', '--masked_only', dest='masked_only', action='store_true', help='Do not write in output file the variables not processed')
		parser.add_argument('-r', '--resolution', choices=cls.resolutions, default=cls.resolution, help='Mask resolution (choices: {%(choices)s}, default: %(default)s)')
		parser.add_argument('-m', '--mode', choices=cls.modes, default=cls.mode, help='Way to decide if a grid point is masked (choices: {%(choices)s}, default: %(default)s)')
		parser.add_argument('-T', '--thresholds', default=cls.thresholds, type=float, nargs=2, help='Thresholds when mode=intersect (default: %(default)s)')
		parser.add_argument('-R', '--reverse', action='store_true', help='Reverse the mask (mask outside)')
		netcdf_group = parser.add_mutually_exclusive_group()
		netcdf_group.add_argument('-3', '--netcdf3', help='Output file with NetCDF3 format (default is NetCDF4)')
		netcdf_group.add_argument('-L', '--netcdf4level', type=int, default=cls.netcdf4level, help='Set NetCDF compression level (0-9, default: %(default)s)')
		return parser
	
	@classmethod
	def main(cls, argv=None):
		parser = cls.get_argparser()
		cls.get_logger().add_argparser_options(parser)
		args = parser.parse_args(argv)
		cls.get_logger().apply_class_argparser_options(args)
		obj = cls(**vars(args))
		obj.run()
		return 0 if obj.error_counter == 0 else 1
	
	def __init__(self, **kwargs):
		
		bases.Object.__init__(self)
		
		# Set attribute from arguments
		# Refer to the main method option parser for details
		for a,v in kwargs.items():
			if v is not None:
				setattr(self, a, v)
		
		# Override error which is setup at runtime in Object.__init__
		self.error_counter = 0
		def error(*args, **kwargs):
			self.error_counter += 1
			return super(self.__class__, self).error(*args, **kwargs)
		self.error = error
	
	def run(self):
	
		if self.netcdf3: io.netcdf3()
		else: io.netcdf4(level=self.netcdf4level)
		cdms2.setAutoBounds(0)

		self.notice('Masking %s to %s', self.inpfile, self.outfile)
		inpfile = cdms2.open(self.inpfile)
		if os.path.isfile(self.outfile):
			if os.path.samefile(self.inpfile, self.outfile):
				raise Exception('Cannot use same input and output file')
			if self.overwrite:
				os.remove(self.outfile)
			else:
				raise Exception('Output file already exists and overwriting is not requested')
		outfile = cdms2.open(self.outfile, 'w')

		# Copy global attributes
		for a,v in inpfile.attributes.items():
			setattr(outfile, a, v)

		# Keep grid masks in memory to improve performances (TODO: use a file cache as for basemaps ?)
		mask_cache = dict()
		stats = True # self.is_verbose()

		# Iterate over input variables
		for varid,filevar in inpfile.variables.items():
			
			# Process only gridded/specified variables
			grid = filevar.getGrid()
			maskit = grid is not None
			if maskit and self.included_variables and varid not in self.included_variables: maskit = False
			if maskit and self.excluded_variables and varid in self.excluded_variables: maskit = False
			if not maskit and self.masked_only:
				self.verbose('Ignoring %s: variable has no grid and masked_only was specified')
				continue
			
			self.logger.info(bases.psinfo())
			self.notice('Processing variable: %s, grid: %s', varid, grid.shape if grid else None)
			# NOTE: With scalar variables, memvar could be numpy.<type> instead of cdms2.tvariable.TransientVariable
			memvar = filevar()
			self.info(bases.describe(memvar, stats=stats))
			self.logger.info(bases.psinfo())
			
			if maskit:
				# Build the mask, check if it is already cached ?
				cache_id = id(grid)
				self.info('Get mask for grid %s', cache_id)
				mask = mask_cache.get(cache_id, None)
				if mask is None:
					self.notice('Loading mask: %s', dict(resolution=self.resolution, mode=self.mode, thresholds=self.thresholds, reverse=self.reverse))
					mask = masking.polygon_mask(grid, self.resolution, mode=self.mode, thresholds=self.thresholds)
					if self.reverse: mask = ~mask
					mask_cache[cache_id] = mask
					self.info(bases.describe(mask, stats=stats))
					self.logger.info(bases.psinfo())
				
				# Mask variable
				# TODO: check/handle dimensions count and order
				self.notice('Masking variable: %s, mask: %s', memvar.shape, mask.shape)
				mask = MV2.resize(mask, filevar.shape)
				memvar[:] = MV2.masked_where(mask, memvar)
				self.info(bases.describe(memvar, stats=stats))
			
			# Special scalar case which could fail if directly written (because fill_value is None)
			if not filevar.shape:
				fill_value = filevar.getMissing()
				if fill_value is None:
					fill_value = -memvar
				memvar = cdms2.createVariable(
					memvar, id=filevar.id, shape=(), typecode=filevar.typecode(),
					fill_value=fill_value, attributes=filevar.attributes)
			
			# Write masked variable to output file
			self.notice('Writing variable to file')
			outfile.write(memvar)

		self.logger.info(bases.psinfo())
		outfile.close()
		inpfile.close()

if __name__ == '__main__':
	sys.exit(Mask().main())
