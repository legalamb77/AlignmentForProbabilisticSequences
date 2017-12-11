'''
@author: Theodore Morley
'''
import numpy as np
import sys

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

    def __init__(self, s, e, pi, t_probs, e_probs):
        self.states = s
        self.emissions = e
        self.initial_probs = pi
        self.transition_probs = t_probs
        self.emission_probs = e_probs
    
    def check_val(self, val):
        if val in self.transition_probs:
            return self.transition_probs[val]
        else:
            return 0

    def viterbi(self, sequence):
        '''
        Takes input in the form of:
            Sequence: An iterable composed of elements from the emissions list.
            #States: A list of states.
            #Emissions: A list of possible emissions.
            #Pi: A dictionary mapping from states to numerical probabilities.
            #Emission_probs: A dictionary mapping from tuples of the form (state, emission) to the probability of said emission inside said state.
            #Transition_probs: A dictionary mapping from tuples of the form (state1, state2) to the probability of the transition from state1 to state2.
        Returns the most likely sequence of hidden states under the assumption of the given model.
        '''
        # Initialize a 2 dimensional structure of size states*emissions. For ease of access we will use a list of dictionaries
        vit = []
        for c in sequence:
            # sanity check
            if c not in self.emissions:
                print("Error! This sequence does not match the emissions given.")
            else:
                vit.append(dict())
        for state in self.states:
            # Fill the first column using the initial probabilities and the emission probabilities
            vit[0][state] = (self.initial_probs[state]*self.emission_probs[(state, sequence[0])], None)
        # Now we can perform the computation of the middle columns in the table
        # We start our range at 1 because we have already filled out the 0th column
        for time in range(1, len(sequence)):
            print(time)
            for state in self.states:
                print(state)
                considered_values = [self.emission_probs[(state, sequence[time])]*self.check_val((state, prev_state))*vit[time-1][prev_state][0] for prev_state in vit[time-1].iterkeys()]
                maximum_arg = max(considered_values)
                backpointer = vit[time-1].keys()[considered_values.index(maximum_arg)]
                vit[time][state] = (maximum_arg, backpointer)
        # Now that the table is filled, we can trace our backpointers back to create our path.
        most_likely_seq = []
        final_max = None
        for state in self.states:
            if vit[len(sequence)-1][state][0] > final_max:
                final_max = state
        max_prob = vit[len(sequence)-1][final_max][0]
        for time in range(len(sequence), 0, -1):
            most_likely_seq.append(final_max)
            final_max = vit[time-1][final_max][1]
        # Final step: reverse the states, and then return.
        most_likely_seq.reverse()
        return (most_likely_seq, max_prob)

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

'''
This function assembles an hmm for the genome sequence which is similar to a simplified profile hmm
It takes into account only insertions and matches, so that we can quickly perform an essentially gapless query check
If we use relatively small word sizes, we can find regions likely to have produced the given query sequence.
The emission probabilities for the match states are set to the probability distribution specified in the genome matrix.
For simplicity, the insertion states have an emission distribution equal to the uniform distribution over the alphabet provided.
'''
def assemble_nodel_hmm(genome_mat, alphabet, tr, window, insertions):
    # tr = [match to match, match to insertion, insertion to self, insertion to match]
    # Needs to take the genome matrix, and construct the following:
    # Also needs to model the irrelevant states -- the ones that we ignore before and after the aligned query
    # States: A list of all possible states
    states = []
    matches = []
    #insertions = []
    for i in range(len(genome_mat[alphabet[0]])):
        matches.append("m"+str(i))
        #if i<len(genome_mat[alphabet[0]]) - 1:
        #    insertions.append("i"+str(i))
    states = matches#+insertions
    # Emissions: A list of all possible emissions (equal to the alphabet)
    emissions = list(alphabet)
    # Initial Probababilities: A uniform distribution over the match states
    pi = dict()
    for state in matches:
        pi[state] = 1.0/float(len(matches))
    #for i_state in insertions:
    #    pi[i_state] = 0
    # Transition probabilities:
    transition_probs = dict()
    for i in range(len(matches)):
        if i<len(matches)-1:
            transition_probs[(matches[i],matches[i+1])] = tr[0]
    #for i in range(len(insertions)):
    #    transition_probs[(matches[i],insertions[i])] = tr[1]
    #    transition_probs[(insertions[i],insertions[i])] = tr[2]
    #    if i<len(matches)-1:
    #        transition_probs[(insertions[i], matches[i+1])] = tr[3]
    # Emission Probabilities: 
    emission_probs = dict()
    for state in matches:
        for char in alphabet:
            emission_probs[(state, char)] = float(genome_mat[char][int(state[1])])
    #for i_state in insertions:
    #    for char in alphabet:
    #        emission_probs[(i_state,char)] = 1.0/len(alphabet)
    return hmm(states, emissions, pi, transition_probs, emission_probs)

def ungapped_eval(start, stop, aligned, score, genome_mat):
    val = 0.0
    pos = start
    for char in aligned:
        val += score[0]*genome_mat[char][pos] + score[1]*(1-genome_mat[char][pos])
        pos += 1
    return val

'''
Performs a heuristic alignment between a query sequence and a genome, based off a given matched zone.
Uses a simulated annealing technique in order to approximate a good alignment, taking the best move at any state and storing that configuration, before
randomly adjusting according to the acceptance probability function
'''
def align(A_seq, B_seq, scoring):
    alignment[['a']['b']]
    return (0, alignment)

'''
Input:
    query: A query sequence
    k: An integer word size
output:
    A list of all k-length words in query
'''
def make_k_list(query, k):
    all_k = []
    pos = 0
    if k>len(query):
        print("Error! k length longer than query!")
        return []
    else:
        all_k.append((query[:k], (0,k-1)))
        pos = k
    while pos+1 <= len(query):
        all_k.append((all_k[-1][0][1:]+query[pos], (pos-k+1, pos)))
        pos+=1
    return all_k

def bumpdown(current_considered, replace_ind, new):
    bumped = []
    for i in range(len(current_considered)):
        if i<replace_ind:
            bumped.append(current_considered[i])
        elif i==replace_ind:
            bumped.append(new)
        else:
            bumped.append(current_considered[i-1])
    return bumped

'''
Input:
    genome_mat: The entire matrix of the genome
    sequence: The entire query sequence
    k and n: The window size to be searched is k -- The query will be broken in to k length words, and the n highest scoring ungapped 
    -alignment according to the eval function will be returned along with their positions to be examined more closely
    score: A tuple of shape (match, mismatch), where match and mismatch are integer scores for those situations
Output: A list of items in the following format
    (score, (genome_start, genome_stop), (aligned_word, (query_start, query_stop)))
    Where the list is arranged from highest to lowest score, and is of length n.
'''
def n_best_ungapped(genome_mat, sequence, n, k, score):
    #considered = []
    print("Finding n best ungapped alignments of length k...")
    nbest = []
    query_klist = make_k_list(sequence, k)
    positions = [(i, i+k) for i in xrange(len(genome_mat[sequence[0]])) if i+k<len(genome_mat[sequence[0]])]
    for pair in positions:
        for qk in query_klist:
            #considered.append((ungapped_eval(pair[0], pair[1], qk[0], score, genome_mat), pair, qk))
            if len(nbest)<n:
                nbest.append((ungapped_eval(pair[0], pair[1], qk[0], score, genome_mat), pair, qk))
            else:
                nbest = sorted(nbest, key=lambda x: x[0], reverse = True)
                potential = (ungapped_eval(pair[0], pair[1], qk[0], score, genome_mat), pair, qk)
                check = 0
                while check<n:
                    if potential[0] > nbest[check][0]:
                        nbest = bumpdown(nbest, check, potential)
                        break
                    else:
                        check+=1
    return nbest#sorted(considered, key=lambda x: x[0])[len(considered)-n:]

def full_search(genome_file, probs_file, query_file, alphabet, initial_n, initial_k, scoring, insertions, tr_probs, hmm_window, hmm_n):
    print("Beginning search...")
    # First, assemble the genome matrix from the provided files
    genome_mat = load_data(genome_file, probs_file, alphabet)
    # Second, read in the query string
    qf = open(query_file, 'r')
    query = qf.read()
    qf.close()
    # Next narrow our search space using a heuristic gapless alignment with initial n and k
    nbest = n_best_ungapped(genome_mat, query, initial_n, initial_k, scoring)
    # Then create simplified HMM's in a predefined window size around the initial promising regions. Insertions are optionally modelled, but deletions are not.
    hmms = []
    for candidate in nbest:
        window = #TODO
        hmms.append((assemble_nodel_hmm(genome_mat, alphabet, tr_probs, window, insertions)), window)
    # Run the viterbi algorithm on each HMM for the query sequence
    candidate_alignments = []
    for hmm, wind in hmms:
        candidate_alignments.append((hmm.viterbi(query),wind))
    # After that, select the best of the hmm outputs (using hmm_n), and build gapped alignments using simulated annealing techniques
    # Finally, report the best of the gapped alignments as our final score
    print("Search complete!")

if __name__ == "__main__":
    print('test')
    alph = ['A','T','G','C']
    t = [0.9, 0.1, 0.7, 0.3]
    genome_matrix = load_data(sys.argv[1],sys.argv[2],alph)
    for v in n_best_ungapped(genome_matrix, 'AAAGGTTCCATCGAAAACCCGGGATCAAACCCTTTTTTGGAGGAGGAGGTAC', 3, 8, (5, -4)):
        print(v)
    #hmm = assemble_nodel_hmm(genome_matrix, alph, t)
    #print(hmm.viterbi('CC'))
    #gf = open(sys.argv[1])
    #text = gf.read()
    #gf.close()
    #make_k_list(text, 8)

