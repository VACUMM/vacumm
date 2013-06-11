from utils import *


class TestSequenceFunctions(VCTestCase):

    for test_name in [
        'test_grid_coord2slice_acad',
#        'test_grid_coord2slice_real',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
