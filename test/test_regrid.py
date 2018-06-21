from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_regrid_curvedinterpolator',
        'test_regrid_extend',
        'test_regrid_fortran_bilin',
        'test_regrid_fortran_bilin2dto1d',
        'test_regrid_fortran_bilin2dto1dc',
        'test_regrid_fortran_cellerr1d',
        'test_regrid_fortran_closest2d',
        'test_regrid_fortran_curv2rect',
        'test_regrid_fortran_curv2rel',
        'test_regrid_fortran_curv2rel_single',
        'test_regrid_fortran_dstwgt',
        'test_regrid_fortran_dstwgt2dto1d',
        'test_regrid_fortran_dstwgt2dto1dc',
        'test_regrid_fortran_extrap1d',
        'test_regrid_fortran_interp1d',
        'test_regrid_fortran_interp1dx',
        'test_regrid_fortran_interp1dxx',
        'test_regrid_fortran_linear4dto1d',
        'test_regrid_fortran_linear4dto1dxx',
        'test_regrid_fortran_nearest2d',
        'test_regrid_fortran_nearest2dto1d',
        'test_regrid_fortran_nearest2dto1dc',
        'test_regrid_grid2xy',
        'test_regrid_griddata',
        'test_regrid_xy2xy',
        'test_regrid_kriging_cloud_split',
        'test_regrid_kriging_krig',
        'test_regrid_kriging_parallel',
        'test_regrid_kriging_regrid',
        'test_regrid_kriging_simple_farvalue',
        'test_regrid_kriging_simple_obserr',
        'test_regrid_kriging_variogram',
        'test_regrid_regrid1d',
        'test_regrid_regrid2d',
        'test_regrid_transect_curvgrid',
        'test_regrid_transect_mld',
        'test_regrid_transect_xyt',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
