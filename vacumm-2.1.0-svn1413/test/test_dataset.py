from utils import *


class TestSequenceFunctions(VCTestCase):

    for test_name in [
        'test_dataset_get_temp_mfs',
        'test_dataset_get_temp_enor',
        'test_dataset_asvar_mfs',
        'test_dataset_asvar_enor',
        'test_dataset_get_dens',
        'test_dataset_plot_transect_menor',
        'test_dataset_plot_transect_mfs',
        'test_dataset_plot_hsection_menor',
        'test_dataset_plot_hsection_mfs',
        'test_dataset_get_uvgbt_menor',
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
