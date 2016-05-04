import re

class fasta:
    ''' read fasta format file'''
    seq = {}
    file = ''
    def __init__(self, file):
        self.file = file
        currentId = ''
        with open(file, 'r') as f:
            for i in f:
                i = i.strip()
                if re.search('^>', i):
                    currentId = i
                    if i not in self.seq.keys():
                        self.seq[i] = ""
                else:
                    self.seq[currentId] = self.seq[currentId] + i
        f.close()

    def getIds(self):
        return self.seq.keys()

    def getSeq(self, id):
        return self.seq[id]


## test
if __name__ == '__main__':
    file = '/Users/shichenwang/Downloads/t.fasta'
    f = fasta(file)

    for id in f.getIds():
        print id
        print f.getSeq(id)