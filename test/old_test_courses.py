from utils import *

method_template = """def test_courses_{0}(self):
    execfile(self.get_path('courses_{0}'))
    self.handle_result(locals().get('result',None))
"""

courses_dir = os.path.realpath(os.path.join(cur_dir,'../scripts/courses'))


class TestSequenceFunctions(VCTestCase):

    def get_path(self, test_name):
        """Get the path to the test script"""
        return os.path.join(courses_dir, test_name+'.py')
    
    for test_name in [
        'python',
        'masking',
        'cdat_bases', 
        'cdat_tools', 
        'interp', 
        'grids', 
        'io_netcdf', 
        'io_other', 
        'numpy', 
        'plots_winds_quiver', 
        'time', 
        'regrid', 
        'bigfiles', 
        'advanced_cf', 
        'advanced_cfgm', 
        'advanced_remote', 
        'advanced_sigma', 
        'advanced_log', 
        'bathy_xyz', 
        'bathy_gridded', 
        'shorelines', 
        'phys_tide', 
        'phys_dyn', 
        ]:
        exec(method_template.format(test_name))
        
        
if __name__ == '__main__':
    unittest.main()
