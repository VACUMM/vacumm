from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_regrid_fortran_extrap1d',
        'test_regrid_fortran_interp1d',
        'test_regrid_fortran_interp1dx',
        'test_regrid_fortran_interp1dxx',
        'test_regrid_fortran_cellerr1d',
        'test_regrid_fortran_bilin',
        'test_regrid_fortran_dstwgt',
        'test_regrid_fortran_nearest2dto1d',
        'test_regrid_fortran_bilin2dto1d',
        'test_regrid_fortran_dstwgt2dto1d',
        'test_regrid_fortran_curv2rect',
        'test_regrid_fortran_nearest2dto1dc',
        'test_regrid_fortran_bilin2dto1dc',
        'test_regrid_fortran_dstwgt2dto1dc',
        'test_regrid_regrid1d',
        'test_regrid_curvedinterpolator',
        'test_regrid_regrid2d',
        'test_regrid_extend',
        'test_regrid_kriging_cloud_split',
        'test_regrid_kriging_variogram',
        'test_regrid_kriging_krig',
        'test_regrid_kriging_parallel',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
