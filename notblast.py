'''
@author: Theodore Morley
'''
import numpy as np
import sys
import math

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
        return np.random.choice(ts, p=ps)
    
    '''
    Returns an emission from the distribution on this state
    '''
    def emit(self):
        es, ps = zip(*self.emissions)
        return np.random.choice(es, p=ps)

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
            #print(time)
            for state in self.states:
                #print(state)
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
def assemble_nodel_hmm(genome_mat, alphabet, tr, insertions):
    print("Building an hmm...")
    # tr = [match to match, match to insertion, insertion to self, insertion to match]
    # Needs to take the genome matrix, and construct the following:
    # Also needs to model the irrelevant states -- the ones that we ignore before and after the aligned query
    # States: A list of all possible states
    states = []
    matches = []
    insertions = []
    for i in range(len(genome_mat[alphabet[0]])):
        matches.append("m"+str(i))
        if insertions:
            if i<len(genome_mat[alphabet[0]]) - 1:
                insertions.append("i"+str(i))
    states = matches+insertions
    # Emissions: A list of all possible emissions (equal to the alphabet)
    emissions = list(alphabet)
    # Initial Probababilities: A uniform distribution over the match states
    pi = dict()
    for state in matches:
        pi[state] = 1.0/float(len(matches))
    if insertions:
        for i_state in insertions:
            pi[i_state] = 0
    # Transition probabilities:
    transition_probs = dict()
    for i in range(len(matches)):
        if i<len(matches)-1:
            transition_probs[(matches[i],matches[i+1])] = tr[0]
    if insertions:
        for i in range(len(insertions)):
            transition_probs[(matches[i],insertions[i])] = tr[1]
            transition_probs[(insertions[i],insertions[i])] = tr[2]
            if i<len(matches)-1:
                transition_probs[(insertions[i], matches[i+1])] = tr[3]
    # Emission Probabilities: 
    emission_probs = dict()
    for state in matches:
        for char in alphabet:
            emission_probs[(state, char)] = float(genome_mat[char][int(state[1])])
    if insertions:
        for i_state in insertions:
            for char in alphabet:
                emission_probs[(i_state,char)] = 1.0/len(alphabet)
    return hmm(states, emissions, pi, transition_probs, emission_probs)

def ungapped_eval(start, stop, aligned, score, genome_mat):
    val = 0.0
    pos = start
    for char in aligned:
        val += score[0]*genome_mat[char][pos] + score[1]*(1-genome_mat[char][pos])
        pos += 1
    return val

'''
score: (match, mismatch, gap_open, gap_general)
'''
def gapped_eval(start, aligned, score, genome_mat):
    val = 0.0
    pos = start
    ingap = False
    for char in aligned:
        if not ingap and char == "_":
            ingap = True
            val += score[2]
        elif ingap and char != "_":
            ingap = False
        if ingap:
            val += score[3]
        else:
            val += score[0]*genome_mat[char][pos] + score[1]*(1-genome_mat[char][pos])
        pos+=1
    return val

def mutate(genome_mat, query, start_g, position, fshift):
    num=1
    if fshift:
        num=3
    for i in range(num):
        query.insert(position, "_")
    start_g = max(0, start_g-num)
    stop = start_g+len(query)
    return (query, start_g, stop)

'''
Performs a heuristic alignment between a query sequence and a genome, based off a given matched zone.
Uses a simulated annealing technique in order to approximate a good alignment, taking the best move at any state and storing that configuration, before
randomly adjusting according to the acceptance probability function
Input:
    genome_mat:
    q_seq:
    scoring: (Match, mismatch, gap_open, gap_general)
    gap_probs: [Probability of len(3) gap, probability of len(1) gap]
'''
def anneal_align(genome_mat, q_seq, scoring, gap_probs, start, T):
    # On each step mutate the alignment by adding either a 3 length gap or a 1 length gap (w/probability defined by gap probs)
    # This is uniformly distributed across positions in the alignment
    # If the resulting score is better, take it, store it, and continue
    # If it is worse, accept it with probability e^(-d/T), where d is (oldscore-newscore), and T is the temperature parameter.
    print("Performing annealing alignment...")
    temperature = True
    q_seq = list(q_seq)
    current_best =final_best= (ungapped_eval(start, start+len(q_seq), q_seq, scoring, genome_mat), (q_seq[:], start, start+len(q_seq)))
    while temperature:
        possible_slots = [1.0/float(len(q_seq)) for x in q_seq]
        pick = np.random.choice(np.arange(start, start+len(q_seq)), p=possible_slots)
        gapThree = np.random.choice([True, False], p=gap_probs)
        candidate = mutate(genome_mat, current_best[1][0], start, pick-start, gapThree)
        candidate_score = gapped_eval(candidate[1], candidate[0], scoring, genome_mat)
        if candidate_score > current_best[0]:
            current_best = (candidate_score, candidate)
        else:
            temperature = math.exp(-(current_best[0]-candidate_score)/float(T))>np.random.random_sample()
            T-=.1
            if temperature:
                current_best = (candidate_score, candidate)
        if current_best[0] > final_best[0]:
            final_best = current_best
    return final_best

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

def getNum(viterbi_result):
    val = viterbi_result[0][0][1:]
    return int(val)

def full_search(genome_file, probs_file, query, res_file, alphabet, initial_n, initial_k, scoring, insertions, tr_probs, hmm_window, hmm_n, gap_probs, temp):
    print("Beginning search...")
    # First, assemble the genome matrix from the provided files
    genome_mat = load_data(genome_file, probs_file, alphabet)
    # Second, read in the query string
    #qf = open(query_file, 'r')
    #query = qf.read().rstrip()
    #qf.close()
    # Next narrow our search space using a heuristic gapless alignment with initial n and k
    nbest = n_best_ungapped(genome_mat, query, initial_n, initial_k, scoring)
    # Then create simplified HMM's in a predefined window size around the initial promising regions. Insertions are optionally modelled, but deletions are not.
    hmms = []
    for candidate in nbest:
        small_genome_mat = dict()
        start = max(candidate[1][0]-hmm_window, 0)
        stop = min(candidate[1][1]+hmm_window, len(genome_mat[alphabet[0]]))
        for char in alphabet:
            small_genome_mat[char] = genome_mat[char][start:stop]
        hmmres = assemble_nodel_hmm(small_genome_mat, alphabet, tr_probs, insertions)
        window = (start, stop)
        hmms.append((hmmres, window))
    # Run the viterbi algorithm on each HMM for the query sequence and score the output
    candidate_alignments = []
    for hmm, wind in hmms:
        vit_res = hmm.viterbi(query)
        ung_score = ungapped_eval(wind[0]+getNum(vit_res), len(vit_res[0][0])+getNum(vit_res), query, scoring, genome_mat)
        candidate_alignments.append((vit_res, (wind[0]+getNum(vit_res), len(vit_res[0][0])+getNum(vit_res)), ung_score))
    # After that, select the best of the hmm outputs (using hmm_n), and build gapped alignments using simulated annealing techniques
    best_hmm = sorted(candidate_alignments, key = lambda x:x[0][1], reverse=True)[:hmm_n]
    # Use this to replace the candidate alignments
    gapped_alignments = []
    for c_hmm in best_hmm:
        gapped_alignments.append(anneal_align(genome_mat, query, scoring, gap_probs, c_hmm[1][0], temp))
    # Finally, report the best of the gapped alignments as our final score, and write the results to the specified file.
    gapped_alignments = sorted(gapped_alignments, key = lambda x: x[0], reverse=True)
    out = open(res_file, 'w')
    for a in gapped_alignments:
        out.write("Result:\n")
        for item in a:
            out.write(str(item)+"\n")
        out.write("\n\n")
    return gapped_alignments[0]
    print("Search complete!")

if __name__ == "__main__":
    alph = ['A','T','G','C']
    t = [0.9, 0.1, 0.7, 0.3]
    gp = [0.99, 0.01]
    qf = open("query.txt", 'r')
    count = 0
    for line in qf:
        query = line.rstrip()
        full_search(sys.argv[1], sys.argv[2], query,"fullres"+str(count)+".txt", alph, 6, 15, (5, -4, -3, -1), True, t, 200, 6, gp, 20)
        count += 1
    qf.close()

