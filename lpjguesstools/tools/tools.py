"""FILE tools.tools.py

Christian Werner, SENCKENBERG Biodiversity and Climate Research Centre (BiK-F)
email: christian.werner@senkenberg.de
2017/11/07
"""

from collections import OrderedDict

def enum(*sequential, **named):
    named.update(dict(unclassified=99))    
    enums = OrderedDict(zip(sequential, range(len(sequential))), **named)
    _enums = enums.copy()
    reverse = dict((value, key) for key, value in enums.items())
    enums['content'] = [reverse[x] for x in _enums.values()]
    enums['items'] = [x for x in _enums.values()]    
    enums['reverse_mapping'] = reverse
    return type('Enum', (), enums)


