from utils import *


class TestArakawa(VCTestCase):

    for test_name in [
        'test_arakawa_interp',
        'test_arakawa_gridtransfer',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
