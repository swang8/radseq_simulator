def reverse(seq):
    rev=[]
    l = len(seq)
    while l > 0:
        l = l - 1
        rev.append(seq[l])
    return ''.join(rev)

def reverse_compliment(seq):
    comp = {'A': "T", 'T':'A', 'G':'C', 'C':'G'}
    revc = []
    l = len(seq)
    while l > 0:
        l = l -1
        s = seq[l]
        if s in comp:
            rev.append(comp[s])
        else:
            rev.append('N')
    return ''.join(revc)

