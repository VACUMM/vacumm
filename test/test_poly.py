from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_poly_create_polygon',
        'test_poly_create_polygons',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
