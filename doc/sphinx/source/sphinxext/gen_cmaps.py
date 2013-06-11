# -*- coding: UTF-8 -*-

# Update the colormaps

import os
import vacumm.misc.color
from sphinx.util.console import bold

def gen_cmaps(app, rstfile):
    
    # Generate specs
    cmaps = {}
    # - direct
    for cmap_name in vacumm.misc.color.cmaps_act():
        cmaps[cmap_name] = cmap_name, 'misc-color-'+cmap_name
    # - manual
    for cmap_funcname, cmap_kwargs, savefigs in [
        ('cmap_jets', dict(stretch=0.6), 'misc-color-vacumm_jets+60'), 
        ('cmap_jets', dict(stretch=-0.6), 'misc-color-vacumm_jets-60'), 
        ('cmap_magic', dict(n=10), 'misc-color-vacumm_magic-n10'), 
        ('cmap_magic', dict(anomaly=True), 'misc-color-vacumm_magic-anom'), 
        ('cmap_magic', dict(positive=True), 'misc-color-vacumm_magic-pos'),
        ('cmap_magic', dict(negative=True), 'misc-color-vacumm_magic-neg'), 
        ]:
        cmap_name = cmap_funcname.replace('cmap_', 'vacumm_')
        cmap_name += '('+','.join([('%s=%s'%item) for item in cmap_kwargs.items()])+')'
        cmaps[cmap_name] = (cmap_funcname, cmap_kwargs), savefigs

    # Plot individual colormaps
    from matplotlib import rc, rcdefaults
    outdir = os.path.join(app.builder.srcdir, os.path.dirname(rstfile))
    for cmap_name in app.builder.status_iterator(
        cmaps.iterkeys(), "generating colormaps... ",
        length=len(cmaps)):
        cmap, savefigs = cmaps[cmap_name]        
        try:
            if not isinstance(cmap, basestring):
                cmap = getattr(vacumm.misc.color, cmap[0])(**cmap[1])
            savefigs = os.path.join(outdir, savefigs)
            vacumm.misc.color.plot_cmap(cmap, savefigs=savefigs, show=False, title=cmap_name, 
                savefigs_verbose=False)
        except:
            pass
            
    # Overview figure
    app.builder.info(bold("generating overview of colormaps..."), nonl=True)
    rc('font', size=8)
    vacumm.misc.color.plot_cmaps(savefigs=os.path.join(outdir, 'misc-color-cmaps'), 
        show=False, figsize=(7, 7), savefigs_verbose=False, aspect=0.1)
    rcdefaults()
    app.builder.info('done')

def check_cmaps(app, env, added, changed, removed):
    """Check if the python file referenced in the rst file has changed
    and regenerate cmaps if so
    """
    rstname = app.config.gen_cmaps_file
    pyfile = None
    for depfile in env.dependencies.get(rstname, ()):
        if depfile.endswith('.py'):
            pyfile = depfile
            break
    else:
        return []
    pyfile = os.path.abspath(os.path.join(env.srcdir, pyfile))
    rstfile = env.doc2path(rstname)
    
    if os.path.getmtime(pyfile)>os.path.getmtime(rstfile):
        gen_cmaps(app, rstfile)
        os.utime(rstfile, None)
        return [rstname]
    return []

def setup(app):
    app.add_config_value('gen_cmaps_file', 'library/misc.color', 'env')
#    app.connect('env-updated', check_cmaps)
    app.connect('env-get-outdated', check_cmaps)
