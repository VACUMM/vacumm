# -*- coding: UTF-8 -*-

# generate a thumbnail gallery of examples
#
# Directly inspired from matplotlib gen_gallery.py
template = """\
{%% extends "layout.html" %%}
{%% set title = "Thumbnail gallery" %%}


{%% block body %%}

<h1>Thumbnail gallery</h1>
<h3>Click to show the code.</h3>
<br/>

<ul>
%s
</ul>
</li>

%s
{%% endblock %%}
"""
link_template = """\
<a href="%s"><img src="%s" border="0" alt="%s"/></a>
"""

header_template = """<div class="section" id="%(subdir)s">\
<h4>%(title)s<a class="headerlink" href="#%(subdir)s" title="Permalink to this headline">¶</a></h4>"""

toc_template = """\
<li><a class="reference internal" href="#%(subdir)s">%(title)s</a></li>"""

import os, glob, re, sys, warnings
import matplotlib.image as image

multiimage_match = re.compile('(.*?)(_\d+)$').match

def make_thumbnail(args):
    image.thumbnail(args[0], args[1], 0.3)

def out_of_date(original, derived):
    return (not os.path.exists(derived) or
            os.stat(derived).st_mtime < os.stat(original).st_mtime)

def gen_gallery(app, doctree):
    if app.builder.name != 'html':
        return

    htmldir = app.builder.outdir
    srcdir = app.builder.srcdir
    
    # From config
    rootdir = app.config.gen_gallery_root
    paths = app.config.gen_gallery_paths
    skips = set(app.config.gen_gallery_skips)
    nmax = app.config.gen_gallery_nmax
    
    # Inits
    thumbnails = {}
    rows = []
    toc_rows = []

    # Loop on subdirs
    for subdir, conf in paths.items():
        
        # Check dirs
        figdir = conf.get('figdir', subdir)
        rstdir = conf.get('rstdir',figdir)
        pngdir = os.path.join(srcdir, figdir)
        thumbdir = os.path.join(htmldir, rootdir, subdir, 'thumbnails')
        if not os.path.exists(thumbdir):
            os.makedirs(thumbdir)

        # Headers
        title = conf.get('title', subdir.title())
        rows.append(header_template % locals())
        toc_rows.append(toc_template % locals())

        # Loop on png files
        data = []
        basenames = []
        for filename in sorted(glob.glob(os.path.join(pngdir, '*.png'))):
            #if filename.endswith("hires.png"):
                #continue

            # Basename
            path, filename = os.path.split(filename)
            basename, ext = os.path.splitext(filename)
            basename = basename.replace('-','.')
            if basename in skips:
                continue
            m = multiimage_match(basename)
            if m is not None:
                basename = m.group(1)
                
            # Check if rst file exists
            if not os.path.exists(os.path.join(srcdir, rstdir, basename+'.rst')):
                continue

            # Check nmax
            if len([True for bn in basenames if bn==basename])>nmax:
                continue
             
            # Is thumbnail out of date?
            orig_path = str(os.path.join(pngdir, filename))
            thumb_path = str(os.path.join(thumbdir, filename))
            if out_of_date(orig_path, thumb_path):
                thumbnails[orig_path] = thumb_path
            
            # Register
            basenames.append(basename)
            data.append((rstdir, basename,
                         os.path.join(rootdir, subdir, 'thumbnails', filename)))


        # Append html rows
        for (rstdir, basename, thumbfile) in data:
            if thumbfile is not None:
                link = '%s/%s.html'%(rstdir, basename)
                rows.append(link_template%(link, thumbfile, basename))

        if len(data) == 0:
            app.warn("no thumbnails were found in %s" % subdir)

        # Close out the <div> opened up at the top of this loop
        rows.append("</div>")

    content = template % ('\n'.join(toc_rows),
                          '\n'.join(rows))

    # Only write out the file if the contents have actually changed.
    # Otherwise, this triggers a full rebuild of the docs
    gallery_path = os.path.join(srcdir, 'templates', 'gallery.html')
    if os.path.exists(gallery_path):
        fh = file(gallery_path, 'r')
        regenerate = fh.read() != content
        fh.close()
    else:
        regenerate = True
    if regenerate:
        fh = file(gallery_path, 'w')
        fh.write(content)
        fh.close()

    # Inform about thumbnail generation
    for key in app.builder.status_iterator(
        thumbnails.iterkeys(), "generating thumbnails... ",
        length=len(thumbnails)):
        image.thumbnail(key, thumbnails[key], 0.3)

def setup(app):
    app.add_config_value('gen_gallery_paths', {}, 'html')
    app.add_config_value('gen_gallery_root', 'gallery', 'html')
    app.add_config_value('gen_gallery_skips', [], 'html')
    app.add_config_value('gen_gallery_nmax', 3, 'html')
    app.connect('env-updated', gen_gallery)

