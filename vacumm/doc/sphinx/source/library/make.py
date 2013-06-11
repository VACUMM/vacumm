"""Creation d'elements necessaire a la doc de l'API"""
print 'Individual colormaps'
from vacumm.misc.color import *
for cmap in cmaps_act():
    print ' -',cmap
    try:
        plot_cmap(cmap,savefigs='misc-color-'+cmap,show=False)
    except:
        pass
plot_cmap(cmap_jets(stretch=0.6), savefigs='misc-color-vacumm_jets+60',title='vacumm_jets(stretch=.6)',show=False)
plot_cmap(cmap_jets(stretch=-0.6), savefigs='misc-color-vacumm_jets-60',title='vacumm_jets(stretch=-.6)',show=False)
plot_cmap(cmap_mg(n=10), savefigs='misc-color-vacumm_magic-n10',title='vacumm_magic(n=10)', show=False)
plot_cmap(cmap_mg(anomaly=True), savefigs='misc-color-vacumm_magic-anom',title='vacumm_magic(anomaly=True)', show=False)
plot_cmap(cmap_mg(positive=True), savefigs='misc-color-vacumm_magic-pos',title='vacumm_magic(positive=True)', show=False)
plot_cmap(cmap_mg(negative=True), savefigs='misc-color-vacumm_magic-neg',title='vacumm_magic(negative=True)', show=False)

print 'All colormaps'
from matplotlib import rc, rcdefaults
rc('font', size=8)
plot_cmaps(savefigs='misc-color-cmaps', show=False, figsize=(7, 7))
rcdefaults()
