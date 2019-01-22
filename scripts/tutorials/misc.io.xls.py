"""Working with excel files with style"""
from builtins import str
from builtins import range
from copy import deepcopy as cp
from xlwt import Workbook, Formula
from vcmqm import xls_style
w = Workbook()

# Define styles
style_header = xls_style(b=1, ha='center', bottom=2, top=2, c=2)
style_left = xls_style(ha='left', i=1)
style_num = xls_style(fmt='0.00')
style_bottom = xls_style( top=2)

# First sheet
# - name
wc = w.add_sheet('Courants')
# - data
data1 = {
    'Four':[1.15646541, 0.25641],
    'Pointe du Raz':[2.56421, 0.1556122],
}
# - width
wc.col(0).width = int(wc.col(0).width * 1.5)
wc.col(1).width = int(wc.col(0).width * 0.7)
# - header
for icol, title in enumerate(('', 'STD [m/s]', 'RMS [%]')):
    wc.write(0, icol, title, style_header)
# - left column and data
for irow, title in enumerate(data1.keys()):
    # title
    wc.write(irow+1, 0, title, style_left)
    # data
    for icol, value in enumerate(data1[title]):
        # RMS in %
        if icol == 1:
            col_style = xls_style(cp(style_num), fmt='0.00%')
        else:
            col_style = style_num
        # write
        wc.write(irow+1, icol+1, str(value), col_style)
# - bottom
for icol in range(3):
    wc.write(irow+2, icol, '', style_bottom)

# Second sheet
# - name
wn = w.add_sheet('Testouille')
# - big title
wn.write_merge(0, 1, 0, 5, u"This is a huge header", xls_style(o=1, s=500))
# - hyperlink
uri = "http://www.ifremer.fr/vacumm"
urn = "VACUMM"
wn.write_merge(3, 3, 1, 10, Formula('HYPERLINK("%s";"%s")'%(uri, urn)),
               xls_style(c=4))
# - levels
wn.write(4, 1, 'level1')
wn.write(5, 1, 'level2')
wn.write(6, 1, 'level2')
wn.write(7, 1, 'level1')
wn.row(4).level = wn.row(7).level = 1
wn.row(5).level = wn.row(6).level = 2
# - dates
from datetime import datetime
wn.write(8, 0, 'Dates:')
wn.write(8, 1, datetime.now(), xls_style(fmt='h:mm:ss AM/PM'))
wn.write(8, 2, datetime.now(), xls_style(fmt='MMM-YY'))

# Write document
w.save('misc.io.xls.xls')
