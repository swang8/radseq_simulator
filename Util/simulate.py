import random

MUT_RATE = 0.001
INDEL_FRAC = 0.15
INDEL_EXTEND = 0.3
MAX_N_RATIO = 0.05

def simulate_error(seq, ERR_RATE = 0.02):
    DNA = ['A', 'T', 'G', 'C']
    err_seq = []
    for i in range(len(seq)):
        if  seq[i] in DNA and random.random() < ERR_RATE:
            ##print 'index:', DNA.index(seq[i])
            position = DNA.index(seq[i])
            subDNA = DNA[0 : position] + DNA[position+1 : ]
            #print subDNA
            err_seq.append(subDNA[random.randint(0, 2)])
        else:
            err_seq.append(seq[i])
    return ''.join(err_seq)


def mutate_genome (id, seq, out_prefix, MUT_RATE = 0.001, INDEL_FRAC = 0.15, INDEL_EXTEND = 0.3, MAX_N_RATIO = 0.05):
    mut_seq = []
    mut_pos = []
    var_file = open(out_prefix + "_var.txt", 'w')
    index = 0
    while index < len(seq):
        nuc = seq[index].upper()
        mut = mutate(nuc)
        #print 'mut:', mut
        if mut == nuc:
            mut_seq.append(nuc)
        else:
            var_file.write(id + '\t' + str(index) + '\t' + nuc + '\t' + mut + '\n')
            mut_pos.append([index, mut])
            if '-' in mut:
                index += len(mut) - 1
            else:
                mut_seq.append(mut)
        index += 1
    var_file.close()
    return [''.join(mut_seq), mut_pos]



def generate_phred_scores(read_length):
    '''
    For the 5 prime part (0 - 2/8); use (30, 10)
    for the middle part (2/8 - 6/8); use (30, 5)
    for the 3 prime part (6/8 -7/8); use (20, 5)
    for the tip of 3 prime part (7/8 -8/8); use (15, 5)
    '''
    range_mu_sd = [[30, 10], [30, 5], [20, 5], [15, 5]]
    phred_scores = []

    for i in range(read_length):
        if i <= (read_length * 2/8):
            (mu, sd) = range_mu_sd[0]
        elif i > (read_length * 2/8) and i <= (read_length * 6/8):
            (mu, sd) = range_mu_sd[1]
        elif i > (read_length * 2/8) and i <= (read_length * 6/8):
            (mu, sd) = range_mu_sd[2]
        else:
            (mu, sd) = range_mu_sd[3]

        score = int(random.normalvariate(mu, sd))
        if (score > 40): score = 40
        if (score < 0): score = 0
        phred_scores.append(score)
    return ''.join([get_chars(e) for e in phred_scores])

def get_chars (p):
    return chr(p + 33)

def mutate(nuc, MUT_RATE = 0.001, INDEL_FRAC = 0.15, INDEL_EXTEND = 0.3, MAX_N_RATIO = 0.05):
    '''
    Given the sequences (seq), simulate the sequencing errors and variations for the sequences.
    MUT_RATE = 0.001
    INDEL_FRAC = 0.15
    INDEL_EXTEND = 0.3
    MAX_N_RATIO = 0.05
    '''
    DNA = ['A', 'T', 'G', 'C', 'N']
    if random.random() < MUT_RATE : ## generate mutation
        mut = ''
        randDNA = DNA[random.randint(0, 4)]
        while randDNA == nuc or randDNA == 'N':
            if  randDNA == "N" and random.random() < MAX_N_RATIO: ## keep the N
                break
            else:
                randDNA = DNA[random.randint(0, 4)]
        ## indel or SNP
        if random.random() < INDEL_FRAC: ## mutation is an indel
            indel_length = random.randint(1, 10) ## max 10bp indel
            indels = []
            if random.random() < 0.5: ## insert
                for i in range(indel_length):
                    indels.append(['A', 'T', 'G', 'C'][random.randint(0,3)])
            else:  # deletion
                for i in range(indel_length):
                    indels.append('-')
            mut = ''.join(indels)
        else:  # SNP
            mut = randDNA
        return mut
    else:  # no mutation happens
        return nuc


if __name__ == '__main__':
    seq = "ATCGATCGTAGCTAGCTAGCTAGTCGATCGTAGCTAGCTA" * 4
    print seq
    K = 0
    while (K < 1000):
        seq_simu = mutate_genome(seq, len(seq))
        K += 1
        if len(seq_simu[1]) == 0:
            continue
        print seq
        for i in seq_simu:
            print i
