from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_grid_coord2slice_acad',
        'test_grid_get_distances',
        'test_grid_subaxis2index',
        'test_grid_resol',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
