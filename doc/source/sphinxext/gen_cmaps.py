# -*- coding: UTF-8 -*-

# Update the colormaps

import os
import vacumm.misc.color
from sphinx.util.console import bold

def gen_cmaps(app, rstfile):

    # Generate specs
    cmaps = {}
    # - direct
    for cmap_name in vacumm.misc.color.cmaps_vacumm():
        cmaps[cmap_name] = cmap_name, app.config.gen_cmaps_prefix+cmap_name
    # - manual
    extra = app.config.gen_cmaps_file

    for cmap_funcname, cmap_kwargs, savefigs in app.config.gen_cmaps_extra_list:
        if not savefigs.startswith(app.config.gen_cmaps_prefix):
            savefigs = app.config.gen_cmaps_prefix+savefigs
        if cmap_kwargs is None: cmap_kwargs = {}
        cmap_name = cmap_funcname.replace('cmap_', 'vacumm_')
        cmap_name += '('+','.join([('%s=%s'%item) for item in cmap_kwargs.items()])+')'
        cmaps[cmap_name] = (cmap_funcname, cmap_kwargs), savefigs

    # Plot individual colormaps
    from matplotlib import rc, rcdefaults
    outdir = os.path.join(app.builder.srcdir, os.path.dirname(rstfile))
    for cmap_name in app.builder.status_iterator(
            cmaps.iterkeys(), "generating colormaps... ",
            length=len(cmaps),
            stringify_func=lambda x: os.path.basename(x),
        ):
        cmap, savefigs = cmaps[cmap_name]
        try:
            if not isinstance(cmap, basestring):
                cmap = getattr(vacumm.misc.color, cmap[0])(**cmap[1])
            savefigs = os.path.join(outdir, savefigs)
            vacumm.misc.color.plot_cmap(cmap, savefigs=savefigs,
                show=False, title=cmap_name, savefigs_verbose=False)
        except:
            pass

    # Overview figure
    app.builder.info(bold("generating overview of colormaps..."), nonl=True)
    rc('font', size=8)
    vacumm.misc.color.plot_cmaps(
        savefigs=os.path.join(outdir, app.config.gen_cmaps_prefix+'cmaps'),
        show=False, figsize=6, savefigs_verbose=False, aspect=0.08, ncol=4,
        close=True)
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
    app.add_config_value('gen_cmaps_prefix', 'misc-color-', 'html')
    app.add_config_value('gen_cmaps_extra_list', [], 'html')
#    app.connect('env-updated', check_cmaps)
    app.connect('env-get-outdated', check_cmaps)
