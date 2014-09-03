from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_plot_curve_o',
        'test_plot_curve_x',
        'test_plot_curve_y',
        'test_plot_curve_t',
        'test_plot_curve_z',
        'test_plot_plot2d_oo',
        'test_plot_plot2d_to',
        'test_plot_section_zo',
        'test_plot_section_quiver',
        'test_plot_hov_tz',
        'test_plot_streamplot',
        'test_plot_core_add_thing',
        'test_plot_taylor', 
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
