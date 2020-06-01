
def suffixArray(s):
    ''' Given T return suffix array SA(T).  Uses "sorted"
        function for simplicity, which is probably very slow.  '''
    satups = sorted([(s[i:], i) for i in range(len(s))])
    ''' extract, return just offsets'''
    return list(map(lambda x: x[1], satups))


def bwtFromSa(t, sa=None):
    ''' Returns BWT(T) for given T by way of the suffix array. '''
    bw = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(t)
    for si in sa:
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si - 1])

    return ''.join(bw), dollarRow  # return string-ized version of list bw
