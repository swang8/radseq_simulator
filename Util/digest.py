from Bio import SeqIO
from Bio.Restriction import *

class Digest:
    fastaFile = ""
    enzymes = []
    fragments = {}

    def __init__(self, enzymes, fastaFile):
        self.enzymes  = enzymes
        self.fastaFile = fastaFile

    def addEnzyme(self, e):
        self.enzymes.append(e)


    def runDigest(self):
        ##self.fragments = {}
        for record in SeqIO.parse(self.fastaFile, "fasta"):
            if record.id not in self.fragments:
                self.fragments[record.id] = {}
            for enz in range(len(self.enzymes)):
                if self.enzymes[enz] in self.fragments[record.id]:
                    continue

                self.fragments[record.id][self.enzymes[enz]] = []
                coords = getattr(Restriction, self.enzymes[enz]).search(record.seq)
                for site in range(len(coords)):
                    print  "%s    %d    %s" % (record.id, coords[site], self.enzymes[enz])
                    self.fragments[record.id][self.enzymes[enz]].append(coords[site])

    def getFragments(self):
        return self.fragments


## for testing
if __name__ == '__main__':
    file = "/Users/shichenwang/Downloads/t.fasta"
    d = Digest(['PstI'], file)
    d.runDigest()
    d.addEnzyme('MluI')
    print 'Add enzyme MluI'
    d.runDigest()
    print d.getFragments()