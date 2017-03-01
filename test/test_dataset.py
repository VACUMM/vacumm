from utils import *


class TSF(VCTestCase):

    for test_name in [
        'test_dataset_get_temp_mfs',
        'test_dataset_get_temp_menor',
        'test_dataset_get_var_seltime',
        'test_dataset_get_variable_levelstring',
        'test_dataset_get_asvar_mfs',
        'test_dataset_get_asvar_menor',
        'test_dataset_get_at',
        'test_dataset_get_depth_menor',
        'test_dataset_get_depth_mfs',
        'test_dataset_get_dens',
        'test_dataset_get_mld',
        'test_dataset_get_uvbt_menor',
        'test_dataset_plot_transect_menor',
        'test_dataset_plot_transect_mfs',
        'test_dataset_plot_hsection_menor',
        'test_dataset_plot_hsection_mfs',
        ]:
        exec(method_template.format(test_name))


if __name__ == '__main__':
    unittest.main()
