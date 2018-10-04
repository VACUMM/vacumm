from utils import *


class TestAxes(VCTestCase):

    for test_name in [
        'test_axes'
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
