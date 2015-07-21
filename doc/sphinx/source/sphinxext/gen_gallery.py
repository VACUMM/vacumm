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
#link_template = """\
#<a href="%s"><img src="%s" border="0" alt="%s"/></a>
#""" # (link, thumbfile, basename)
link_template = """\
<figure>
    <a href="%(link)s"><img src="%(thumbfile)s" border="0" alt="%(basename)s"/></a><br/>
    <figcaption><a href="%(link)s">%(basename)s</a></figcaption>
</figure>
"""

header_template = u"""<div class="section" id="%(subdir)s">\
<h3>%(title)s<a class="headerlink" href="#%(subdir)s" title="Permalink to this headline">¶</a></h3>"""

toc_template = """\
<li><a class="reference internal" href="#%(subdir)s">%(title)s</a></li>"""

import os, glob, re, sys, warnings
import matplotlib.image as image
from matplotlib import rcParams
from codecs import open

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

                # Ok loof for a rst file whose basename matches end for png basename
                rstfiles = filter(
                    lambda fname: basename.endswith(
                        os.path.splitext(os.path.basename(fname))[0]),
                    glob.glob(os.path.join(srcdir, rstdir, '*.rst')))
                if not rstfiles:
                    continue
                basename = os.path.splitext(os.path.basename(rstfiles[0]))[0]


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
                rows.append(link_template%locals())

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
        fh = open(gallery_path, 'r', 'utf8')
        regenerate = fh.read() != content
        fh.close()
    else:
        regenerate = True
    if regenerate:
        fh = open(gallery_path, 'w', 'utf8')
        fh.write(content)
        fh.close()

    # Inform about thumbnail generation
    old = rcParams['image.origin']
    rcParams['image.origin'] = 'upper'
    for key in app.builder.status_iterator(
            thumbnails.iterkeys(), "generating thumbnails... ",
            length=len(thumbnails)):
        image.thumbnail(key, thumbnails[key], 0.3)
    rcParams['image.origin'] = old

def setup(app):
    app.add_config_value('gen_gallery_paths', {}, 'html')
    app.add_config_value('gen_gallery_root', 'gallery', 'html')
    app.add_config_value('gen_gallery_skips', [], 'html')
    app.add_config_value('gen_gallery_nmax', 3, 'html')
    app.connect('env-updated', gen_gallery)

