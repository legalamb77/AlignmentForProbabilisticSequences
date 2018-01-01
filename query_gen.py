import sys
import numpy as np
alph = ['A', 'C', 'G', 'T']

'''
Takes as input the name of the file containing the genome, and the name of the file specifying the probabilities for said genome at each position.
Returns a dict of lists in the following format:
    For each character in the alphabet (Currently assumed to be ACGT), there is a list of length n, where n is the length of the genome
    Each position in this list gives the probability that this position has the given character
'''
def load_data(genome_file, probs_file, alphabet):
    matrix = dict()
    gf = open(genome_file)
    genome = gf.read()
    gf.close()
    pf = open(probs_file)
    probs = pf.read()
    pf.close()
    probs = probs.split()
    genome = list(genome)
    for char in alphabet:
        matrix[char] = []
    for i in range(len(genome)):
        if genome[i] not in alphabet:
            pass
        else:
            others = [c for c in alphabet if c!=genome[i]]
            otherProb = (1-float(probs[i]))/float(len(alphabet)-1)
            matrix[genome[i]].append(float(probs[i]))
            for nucleotide in others:
                matrix[nucleotide].append(otherProb)
    return matrix

def insertion(prob):
    if np.random.random_sample() < prob:
        return True
    else:
        return False

def deletion(prob):
    if np.random.random_sample() < prob:
        return False
    else:
        return True

def substitution(el, prob):
    if np.random.random_sample() < prob:
        return np.random.choice(alph)
    else:
        return el

def main():
    genome_file = sys.argv[1]
    prob_file = sys.argv[2]
    length = int(sys.argv[3])
    sub = float(sys.argv[4])
    insert = float(sys.argv[5])
    delete = float(sys.argv[6])
    num_queries = int(sys.argv[7])
    g_mat = load_data(genome_file, prob_file, alph)
    positions = np.arange(len(g_mat[alph[0]]))
    starts = []
    for i in range(num_queries):
        starts.append(np.random.choice(positions))
    queries = []
    for i in range(num_queries):
        query = []
        current = starts[i]
        for i in range(length):
            pro = []
            for char in alph:
                pro.append(g_mat[char][current])
            query.append(np.random.choice(alph,p=pro))
            current+=1
        queries.append(query)
    #Write unmutated
    out3 = open("unmutated.txt", "w")
    for i in range(len(queries)):
        for char in queries[i]:
            out3.write(char)
        out3.write(", "+str(starts[i]))
        out3.write("\n")
    #insertions
    for i in range(len(queries)):
        for j in range(len(queries[i])):
            if insertion(insert):
                queries[i].insert(j, np.random.choice(alph))
        #deletions
        queries[i] = [x for x in queries[i] if deletion(delete)]
        #substitutions
        queries[i] = [substitution(x, sub) for x in queries[i]]
    out = open("query.txt", "w")
    for i in range(len(queries)):
        for q in queries[i]:
            for char in q:
                out.write(char)
        out.write("\n")
    out.close()
    #out2 = open("params.txt", "w")
    #for i in range(num_queries):
    #    out2.write(str(starts[i])+"\n")
    #out2.close()
    print(len(g_mat[alph[0]]))

if __name__ == "__main__":
    main()
