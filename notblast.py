'''
@author: Theodore Morley
'''
import numpy as np

class state:
    "A state for use in an HMM"

    '''
    The transition probs are a list of tuples in the format (state, probability of transition)
    The same format applies for emission probs, (emission, probability of emission)
    Both of these lists form probability distributions
    '''
    def __init__(self, transition_probs, emission_probs):
        self.transitions = transition_probs
        self.emissions = emission_probs

    '''
    Selects the next state to move to and returns it
    '''
    def step(self):
        ts, ps = zip(*self.transitions)
        return numpy.random.choice(ts, p=ps)
    
    '''
    Returns an emission from the distribution on this state
    '''
    def emit(self):
        es, ps = zip(*self.emissions)
        return numpy.random.choice(es, p=ps)

class hmm:

    '''
    In this case, the probability of emissions and transitions are defined within each state
    '''
    def __init__(self, s, e, pi):
        self.states = s
        self.emissions = e
        self.initial_probs = pi

    '''
    Returns the most likely sequence of states, given this sequence of observations
    '''
    def viterbi(self, observed):
        # for state in
        return True

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
    probs = list(probs)
    genome = list(genome)
    for char in alphabet:
        matrix[char] = []
    for i in range(len(genome)):
        if genome[i] not in alphabet:
            pass
        else:
            others = [c for c in alphabet if c!=genome[i]]
            otherProb = (1-probs[i])/float(len(alphabet)-1)
            matrix[genome[i]].append(probs[i])
            for nucleotide in others:
                matrix[nucleotide].append(otherProb)
    return matrix

'''
'''
def align(A_seq, B_seq, scoring):
    alignment[['a']['b']]
    return (0, alignment)

'''
Input:
    genome_matrix: Result of load_data
    start: The start position of the window being compared to the query, inclusive
    stop: The stop position of the window being compared to the query, inclusive
    query: A list of characters of length stop-start+1.
    n: The number of assignments to the genome_matrix to return
    alphabet: The characters used in the sequence
Output: The n best assignments of characters to the positions in the genome_matrix segment with respect to the query
    For explanation, we want to assume that this region is going to generate the given sequence, so we find the most likely assignments given that assumption.
    We return a list of tuples of the form (assignment, probability)
'''
def n_best_matches(genome_matrix, start, stop, query, n, alphabet):
    return True

'''
Input:
    genome_matrix: The result of the load_data function above
    query: A list of characters in the same alphabet as the genome_matrix
    score_threshold:
    n_threshold:
Output:
    #
'''
def notblast(genome_matrix, query, score_threshold, n_threshold):
    return True


def blast(query, genome, k):
    #Remove sequence repeats
    #Make a k-letter word list of the query sequence
    #List the possible matching words
    #Organize high-scoring words into a search tree
    #Repeat last two steps for all k letter words in query
    #Scan the database sequences for exact matches with the remaining high-scoring wordss
    #Extend the exact matches to high-scoring segment pair
    #List all HSP's in the database higher than the cutoff score
    #Check significance of HSP score
    return True

def prob_align(query, genome):
    #
