## restriction enzymes
from Bio.Restriction import *

class enzymeDict:

    def __init__(self):
        self.re_dict = Restriction_Dictionary.rest_dict

    def ids(self):
        return self.re_dict.keys()

    def getProperty(self, id):
        return self.re_dict[id]

    def checkEnz(self, enz): ## check if the enzyme is in the dictionary.
        return enz in self.re_dict

if __name__ == '__main__':
    ezd = enzymeDict()
    for i in enzymeDict().ids():
        print i
        print ezd.getProperty(i)




