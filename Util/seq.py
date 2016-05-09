def reverse(seq):
    rev=[]
    l = len(seq)
    while l > 0:
        l = l - 1
        rev.append(seq[l])
    return ''.join(rev)

def reverse_complement(seq):
    comp = {'A': "T", 'T':'A', 'G':'C', 'C':'G'}
    revc = []
    l = len(seq)
    while l > 0:
        l = l -1
        s = seq[l]
        if s in comp:
            revc.append(comp[s])
        else:
            revc.append('N')
    return ''.join(revc)

