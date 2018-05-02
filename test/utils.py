import unittest
import os
from matplotlib import use, rcdefaults
use('Agg')
from matplotlib import pyplot as P
#from vcmq import *

rcdefaults()

cur_dir = os.path.abspath(os.path.dirname(__file__))
test_dir = os.path.realpath(os.path.join(cur_dir,'../scripts/test'))

method_template = """def {0}(self):
    P.close()
    execfile(self.get_path('{0}'))
    self.save_figs('{0}')
    self.handle_result(locals().get('result',None))
    P.rcdefaults()
"""

class VCTestCase(unittest.TestCase):

    def get_path(self, test_name):
        """Get the path to the test script"""
        return os.path.join(test_dir, test_name+'.py')

    def save_figs(self, test_name):
        """Save all pending figure to numbered files"""
        pyfile = self.get_path(test_name)
        for i, fig_num in enumerate(P.get_fignums()):
            fig = P.figure(fig_num)
            fig_path = pyfile[:-3] + '_{:02d}.png'.format(i)
            fig.savefig(fig_path)
        P.close()

    def handle_result(self, result=None):
        """Handle result from test scripts"""
        if isinstance(result, dict): result = list(result.items())
        if not isinstance(result, (list, tuple)): return

        # Loop on content
        for key, values in result:
            self.check_single_result(key, values)


    def check_single_result(self, key, values):

#        # Files exist
#        if key=='files':
#            files = values
#            if isinstance(files, basestring):
#                files = [files]
#            for fn in files:
#                self.assertTrue(os.path.exists(fn))
#            return

        # Function that returns True or None if succeeded
        if callable(key):
            if not isinstance(values, (list,tuple)):
                values = [values]
            res = key(*values)
            if key.__module__!='numpy.testing.utils':
                self.assertTrue(res)
            return

        # Assertions
        if key.startswith('assert'):
            if not isinstance(values, (list,tuple)):
                values = [values]
            getattr(self, key)(*values)
            return

