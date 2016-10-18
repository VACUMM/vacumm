from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_color_anamorph',
        'test_color_whiten',
        'test_color_darken',
        'test_color_saturate',
        'test_color_desaturate',
        'test_color_pastelise',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
