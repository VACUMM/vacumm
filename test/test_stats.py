from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_stats_stataccum_single',
        'test_stats_stataccum_dual',
        'test_stats_stataccum_dumpload',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
