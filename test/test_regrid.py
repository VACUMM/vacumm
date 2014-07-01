from utils import *


class TestSequenceFunctions(VCTestCase):

    for test_name in [
        'test_regrid_extend',
        'test_regrid_kriging_cloud_split',
        'test_regrid_kriging_variogram',
        'test_regrid_kriging_krig',
        'test_regrid_kriging_parallel',
        'test_regrid_regrid1d',
        'test_regrid_fortran_extrap1d',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
