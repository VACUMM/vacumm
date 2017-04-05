from vcmq import create_axis, isaxis, coord2slice, cdms2, N


def subaxis2slice(cdaxis, values):
    if cdms2.isVariable(cdaxis):
        return [subaxis2index(cdax, vv) for cdax, vv in
            zip(cdaxis.getAxisList(), values)]
    cdaxis = create_axis(cdaxis)
    ijk = cdaxis.mapIntervalExt((values[0], values[-1], 'cc'))
    if ijk is None:
        return
    return slice(*ijk)


cdaxis = create_axis(N.linspace(0, 11., 17))
subaxis = cdaxis[2:7]

res = subaxis2slice(cdaxis, subaxis)


