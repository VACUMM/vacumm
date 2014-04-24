from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_cdat_regrid_algo',
        'test_cdat_regrid_conserv',
        'test_cdat_regrid_rect2rect',
        'test_cdat_regrid_curv2rect',
        'test_cdat_regrid_regrid2',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
