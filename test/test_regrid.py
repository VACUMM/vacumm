from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_regrid_fortran_extrap1d',
        'test_regrid_fortran_interp1d',
        'test_regrid_fortran_interp1dx',
        'test_regrid_fortran_interp1dxx',
        'test_regrid_fortran_curv2rect',
        'test_regrid_regrid1d',
        'test_regrid_curvedinterpolator',
        'test_regrid_extend',
        'test_regrid_kriging_cloud_split',
        'test_regrid_kriging_variogram',
        'test_regrid_kriging_krig',
        'test_regrid_kriging_parallel',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
