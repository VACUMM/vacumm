"""Test :class:`~vacumm.misc.misc.dictree_get`"""
from vcmq import dicttree_get


dd = {
    'key1':{
        'key11':'val11',
        },
    'key2':'val2',
}

assert dicttree_get(dd, 'key1') == {'key11': 'val11'}
assert dicttree_get(dd, 'key2') == 'val2'
assert dicttree_get(dd, 'key3') == None
assert dicttree_get(dd, 'key1', 'key11') == 'val11'
assert dicttree_get(dd, 'key1', 'key11', 'key111') == 'val11'
