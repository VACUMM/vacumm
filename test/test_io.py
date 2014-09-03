from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_io_shapes',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
