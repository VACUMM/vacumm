# -*- coding: utf8 -*-
"""
Generic tools using markup for html report generation

"""
# Copyright or © or Copr. Ifremer (contributor(s) : Guillaume Charria) (2010)
# 
# guillaume.charria@ifremer.fr
# 
# 
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
# 

from vacumm.markup import markup

__all__ = ['simplereporthtml','simplecomparisonhtml','simplehtml']

def simplehtml(title = None,  header = None,  footer = None,  images_results = None, file_out = None,  intro = None):
    """Quick generation of html report

    simplereporthtml:
    
    - *title*: Title of the html page [default: '--- VACUMM ...  ---']
    - *header*: Text at the beginning of the webpage [default: '____-]
    - *intro*: Introduction text
    - *footer*: Text at the end of the webpage [default: '____']
    - *images_results*: List of results figures to display
    - *file_out*: Name of the html file written [default: 'report_html.html']
    """
    if title is None:
        title = '--- VACUMM ... ---'
    if intro is None:
        intro = ''
    if header is None:
        header = 65*'_'
    if footer is None:
        footer = 65*'_'
    if file_out is None:
        file_out = 'report_html.html'
    #if img is None:

    page = markup.page( )
    page.init( title=title, header=header, footer=footer )
    page.br( )

    # ----
    page.h1('Evaluation de la simulation')
    page.br( )
    page.h2(intro)
    page.div( style="margin:0px;margin-bottom:0px")
    page.img( src=images_results)
    page.div.close()

   # ---- print page    
    f = open(file_out, 'w')
    f.write(str(page))
    f.close()



    return


def simplereporthtml(title = None,  header = None,  footer = None,  images_results = None,  images_control = None,  mimamodel = None,  mimaobs = None,  file_out = None,  intro = None):
    """Quick generation of html report

    simplereporthtml:
    
    - *title*: Title of the html page [default: '--- VACUMM ...  ---']
    - *header*: Text at the beginning of the webpage [default: '____-]
    - *intro*: Introduction text
    - *footer*: Text at the end of the webpage [default: '____']
    - *images_results*: List of results figures to display
    - *images_control*: List of control figures to display
    - *mima*: Text for min and max    
    - *file_out*: Name of the html file written [default: 'report_html.html']
    """


    if title is None:
        title = '--- VACUMM ... ---'
    if intro is None:
        intro = ''
    if header is None:
        header = 65*'_'
    if footer is None:
        footer = 65*'_'
    if file_out is None:
        file_out = 'report_html.html'
    #if img is None:
    
    page = markup.page( )
    page.init( title=title, header=header, footer=footer )
    page.br( )
        
    # ----
    page.h1('Evaluation de la simulation')
    page.br( )
    page.h2(intro)
    page.div( style="margin:0px;margin-bottom:0px")
    page.h4(mimamodel)
    page.h4(mimaobs)
    #page.img( width=600, height=600, src=images_results)  
    page.img(src=images_results) 
    page.div.close()
    page.add('<hr>')
    page.h1('Diagnostics de controle')
    page.div( style="margin:0px;margin-bottom:0px")
    #page.img( width=600, height=600, src=images_control)  
    page.img(src=images_control)
    page.div.close()
    
    
    # ---- print page    
    f = open(file_out, 'w')
    f.write(str(page))
    f.close()
    
    
    
    return


def simplecomparisonhtml(title = None,  header = None,  footer = None, images_cov = None,  images_mean = None,  images_std = None,  images_biasmap = None, images_bias1d = None, file_out = None,  intro = None):
    """Quick generation of html report

    simplecomparisonhtml:
    
    - *title*: Title of the html page [default: '--- VACUMM ...  ---']
    - *header*: Text at the beginning of the webpage [default: '____-]
    - *intro*: Introduction text
    - *footer*: Text at the end of the webpage [default: '____']
    - *images_cov*: List of results figures to display

    - *images_mean*: List of results figures to display

    - *images_std*: List of results figures to display

    - *images_biasmap*: List of results figures to display

    - *images_bias1d*: List of results figures to display

    - *file_out*: Name of the html file written [default: 'report_html.html']
    """


    if title is None:
        title = '--- VACUMM ... ---'
    if intro is None:
        intro = ''
    if header is None:
        header = 65*'_'
    if footer is None:
        footer = 65*'_'
    if file_out is None:
        file_out = 'report_html.html'
    #if img is None:

    page = markup.page( )
    page.init( title=title, header=header, footer=footer )
    page.br( )

    # ----
    page.h1('Evaluation de la simulation')
    page.br( )
    page.h2(intro)
    page.div( style="margin:0px;margin-bottom:0px")
    page.img( width= 600, height=600, src=images_cov)
    page.div.close()
    page.div( style="margin:0px;margin-bottom:0px")
    page.img( width= 500, height=500, src=images_mean)
    page.div.close()
    page.div( style="margin:0px;margin-bottom:0px")
    page.img( width=500, height=500, src=images_std)
    page.div.close()
    page.div( style="margin:0px;margin-bottom:0px")
    page.img( height=600, width=600, src=images_biasmap)
    page.div.close()
    page.div( style="margin:0px;margin-bottom:0px")
    page.img( height=600, width=600, src=images_bias1d)
    page.div.close()
    page.add('<hr>')
    page.div.close()


    # ---- print page    
    f = open(file_out, 'w')
    f.write(str(page))
    f.close()

    return


