from utils import *


class TSMisc(VCTestCase):

    for test_name in [
        'test_misc_casechecker',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
