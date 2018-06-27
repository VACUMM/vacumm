import os
import glob
import subprocess
from sphinx.util.console import bold
from sphinx.util import status_iterator
from vacumm.misc.config import opt2rst


def write_helps(app):

    srcdir = app.env.srcdir

    scripts_patterns = app.config.gen_binhelps_scripts_patterns
    scripts_bindir = os.path.abspath(app.config.gen_binhelps_scripts_bindir)
    rstfilename_format = app.config.gen_binhelps_rstfilename_format
    rstfiles_toupdate = app.config.gen_binhelps_rstfiles_toupdate

    # Guess the script list
    if isinstance(scripts_patterns, str):
        scripts_patterns = [scripts_patterns]
    scripts = []
    for ss in scripts_patterns:
        for script in glob.glob(os.path.abspath(
                os.path.join(scripts_bindir, ss))):

            # Is it executable?
            if not os.path.isfile(script) or not os.access(script, os.X_OK):
                continue

            # Rst file name
            rstfile = os.path.join(srcdir, rstfilename_format.format(
                basename=os.path.basename(script[:-3])))

            # Checks time
            if (os.path.exists(rstfile) and
                    os.stat(script).st_mtime < os.stat(rstfile).st_mtime):
                continue

            scripts.append((rstfile, script))

    for rstfile, script in status_iterator(
            scripts, bold('generating help for bin scripts... '),
            length=len(scripts),
            stringify_func=lambda x: os.path.basename(x[0])):

        # Generate help
        try:
            std,err = subprocess.Popen([script, "-h"], stdout=subprocess.PIPE,
                stderr=subprocess.PIPE).communicate()
        except:
            app.warn('error with'+script)
            continue
        std = std.decode('utf-8', 'replace')
        err = err.decode('utf-8', 'replace')
        if not err is not None: continue

        # Write to file
        dirname = os.path.dirname(rstfile)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        f = open(rstfile, 'w')
        f.write(opt2rst(std))
        f.close()

    # Force documents update with touch
    if scripts:
        if isinstance(rstfiles_toupdate, str):
            rstfiles_toupdate = [rstfiles_toupdate]
        for rstfile in rstfiles_toupdate:
            rstfile = os.path.join(srcdir, rstfile)
            if os.path.exists(rstfile):
                touch(rstfile)

def touch(fname, times=None):
    fhandle = open(fname, 'a')
    try:
        os.utime(fname, times)
    finally:
        fhandle.close()


def setup(app):

    app.add_config_value('gen_binhelps_rstfiles_toupdate', [], '')
    app.add_config_value('gen_binhelps_rstfilename_format',
        'bin/{basename}.help.txt', '')
    app.add_config_value('gen_binhelps_scripts_patterns', ['*.py'], '')
    app.add_config_value('gen_binhelps_scripts_bindir',
        os.path.abspath('../../bin'), '')

    app.connect('builder-inited', write_helps)

    return {'version': '0.1'}

