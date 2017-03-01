from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_atime_interp_clim',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
