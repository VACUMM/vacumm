from utils import *


class TestSequenceFunctions(VCTestCase):

    for test_name in [
        'test_regrid_extend',
        'test_regrid_kriging_coud_split',
        'test_regrid_kriging_variogram',
        'test_regrid_kriging_krig',
        'test_regrid_kriging_parallel',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
