#!/usr/bin/env python

"""EM algorithm for learning rhyming words and rhyme schemes with independent stanzas.
Sravana Reddy (sravana@cs.uchicago.edu), 2011.
"""

import sys, string, random, pickle, os, re
import numpy, math
from collections import defaultdict
import codecs

def load_stanzas(filename):
    """Load raw stanzas from gold standard file"""
    f=codecs.open(filename, 'r', 'utf8').readlines()
    stanzas=[]
    for i, line in enumerate(f):
        line=line.split()
        if i%4==0:
            stanzas.append(line[1:])
    return stanzas

def load_schemes(schemefile):
    """load rhyme schemes from pickled file"""
    schemes=pickle.load(open(schemefile))
    for i in schemes:
        schemes[i]=map(lambda x:map(int, x[0].split()), schemes[i])  #remove freq and convert to list of integers
    return schemes

def get_wordset(stanzas):
    """get all words"""
    words=sorted(list(set(reduce(lambda x, y: x+y, stanzas))))
    return words

def get_rhymelists(stanza, scheme):
    """transform stanza into ordered lists of rhymesets as given by rhyme scheme"""
    rhymelists=defaultdict(list)
    for stanzaword, schemeword in zip(stanza, scheme):
        rhymelists[schemeword].append(stanzaword)
    return rhymelists.values()

def init_uniform_ttable(words):
    """initialize (normalized) theta uniformly"""
    n=len(words)
    t_table=numpy.zeros((n, n+1))  
    uni_prob=1/float(n)
    for c in range(n+1):
        for r in range(n):
            t_table[r, c]=uni_prob   
    return t_table
    
def basic_word_sim(word1, word2):
    """Simple measure of similarity: num of letters in common/max length"""
    common=0.0
    if word1==word2:
        return 1.0
    for c in word1:
        if c in word2:
            common+=1
    return common/max(len(word1), len(word2))

def init_basicortho_ttable(words):
    """initialize probs according to simple measure of orthographic similarity"""
    n=len(words)
    t_table=numpy.zeros((n, n+1)) 

    #initialize P(c|r) accordingly
    for r, w in enumerate(words):
        for c, v in enumerate(words):
            if c<r:
                t_table[r, c]=t_table[c, r]  #similarity is symmetric
            else:
                t_table[r, c]=basic_word_sim(w, v)+0.001  #for backoff
        t_table[r, n]=random.random()  #no estimate for P(r|no history)

    #normalize
    for c in range(n+1):
        tot=float(sum(t_table[:, c]))
        for r in range(n):
            t_table[r, c]=t_table[r, c]/tot

    return t_table

vowels=re.compile('[iye|aou#$3(IYE/A{&QO}VU@!)*<cq0~^KLM123456789WBX]')
celexdir='../../data/celex/CELEX_V2/'  #change to the location of your CELEX directory
epwfile=celexdir+'/english/epw/epw.cd'

def read_celex():
    spam=map(lambda x:x.strip().split('\\'), open(epwfile).readlines())
    spam=map(lambda x:(x[1], x[6].replace('-', '').replace('"', "'")), spam)
    d=defaultdict(list)
    for (word, pron) in spam:
        if "'" in pron:   #can only test words with at least on stressed syllable
            d[word].append(pron)
    return d

def isRhyme(d, w1, w2):
    """check if words rhyme"""
    for p1 in d[w1]:
        #extract only "rhyming portion"
        p1=p1.split("'")[-1]
        m=vowels.search(p1)
        if not m:
            print p1
        p1=p1[m.start():]
        for p2 in d[w2]:
            p2=p2.split("'")[-1]
            m=vowels.search(p2)
            if not m:
                print w2, p2
            p2=p2[m.start():]
            if p1==p2:
                return True
    return False

def init_perfect_ttable(words):
    """initialize (normalized) theta according to whether words rhyme"""
    d=read_celex()
    
    not_in_dict=0
    
    n=len(words)
    t_table=numpy.zeros((n, n+1))  
    
    #initialize P(c|r) accordingly
    for r, w in enumerate(words):
        if w not in d:
            not_in_dict+=1
        for c, v in enumerate(words):
            if c<r:
                t_table[r, c]=t_table[c, r]
            elif w in d and v in d:
                t_table[r, c]=int(isRhyme(d, w, v))+0.001  #for backoff
            else:
                t_table[r, c]=random.random()
        t_table[r, n]=random.random()  #no estimate for P(r|no history)

    print not_in_dict, "of", n, " words are not in CELEX"
    
    #normalize
    for c in range(n+1):
        tot=float(sum(t_table[:, c]))
        for r in range(n):
            t_table[r, c]=t_table[r, c]/tot

    return t_table

def post_prob_scheme(t_table, words, stanza, myscheme):
    """posterior prob of a scheme for a stanza, with prob of every word in rhymelist rhyming with all one before it"""
    myprob=1.0
    n=len(words)
    rhymelists=get_rhymelists(stanza, myscheme)
    plen=len(stanza)
    for rhymelist in rhymelists:
        for i, w in enumerate(rhymelist):            
            r=words.index(w)
            if i==0:  #first word, use P(w|x)
                myprob=myprob*t_table[r, n]
            else:
                for v in rhymelist[:i]:  #history
                    c=words.index(v)
                    myprob*=t_table[r, c]
    if myprob==0 and len(stanza)>30: #probably underflow
        myprob=1e-300
    return myprob

def e_unnorm_post(t_table, words, stanzas, schemes, rprobs):
    """compute posterior prob of rhymescheme for each stanza (expectation step)"""
    probs=[]
    numstanzas=len(stanzas)
    for i, stanza in enumerate(stanzas):
        if i==numstanzas/2:
           print i
        elif i%10==0:
            sys.stdout.write('.')
        stanzaprobs=[]
        myschemes=schemes[len(stanza)]
        for myscheme in myschemes:
            stanzaprobs.append(rprobs[tuple(myscheme)]*post_prob_scheme(t_table, words, stanza, myscheme)) 
        probs.append(stanzaprobs)
    print
    return probs

def e_norm_post(probs):
    """normalize posterior probs"""
    normprobs=[]
    for stanzaprobs in probs:
        tot=sum(stanzaprobs)
        if tot>0:
            normstanzaprobs=map(lambda myprob: myprob/tot, stanzaprobs)
        else:
            normstanzaprobs=stanzaprobs[:]
        normprobs.append(normstanzaprobs)
    return normprobs

def m_frac_counts(words, stanzas, schemes, normprobs):
    """find fractional pseudocounts (maximization step)"""
    n=len(words)
    tc_table=numpy.zeros((n, n+1))
    rprobs=defaultdict(float)
    for stanza, stanzaprobs in zip(stanzas, normprobs):
        myschemes=schemes[len(stanza)]
        for myscheme, myprob in zip(myschemes, stanzaprobs):

            rprobs[tuple(myscheme)]+=myprob  

            rhymelists=get_rhymelists(stanza, myscheme)
            for rhymelist in rhymelists:
                for i, w in enumerate(rhymelist):
                    r=words.index(w)
                    tc_table[r, n]+=myprob
                    for v in rhymelist[:i]+rhymelist[i+1:]:
                        c=words.index(v)
                        tc_table[r, c]+=myprob


    return [tc_table, rprobs]

def m_norm_frac(tc_table, n, rprobs):
    """normalize counts to get conditional probs"""
    t_table=numpy.zeros((n, n+1))

    for c in range(n+1):
        tot=sum(tc_table[:, c])
        if tot==0:
            continue
        for r in range(n):
            t_table[r, c]=tc_table[r, c]/tot

    
    totrprob=sum(rprobs.values())
    for scheme in rprobs:
        rprobs[scheme]=rprobs[scheme]/totrprob
    
    
    return [t_table, rprobs]

def iterate(t_table, words, stanzas, schemes, rprobs, maxsteps):
    """iterate steps 2-5 until convergence, return final t_table"""
    numstanzas=float(len(stanzas))
    data_prob = -10**10
    epsilon=.1
    for ctr in range(maxsteps):
        old_data_prob=data_prob
        
        #E-step
        probs=e_unnorm_post(t_table, words, stanzas, schemes, rprobs)
        
        #estimate total probability
        allschemeprobs=map(sum, probs)

        if 0.0 in allschemeprobs:   #this may happen for very large data on large stanzas, small hack to prevent
            underflows=filter(lambda x:x[2]==0.0, zip(range(len(stanzas)), stanzas, allschemeprobs))
            for underflow in underflows:
                if len(probs[underflow[0]])==1:
                    probs[underflow[0]][0]=1e-300
                    allschemeprobs[underflow[0]]=1e-300
                    print "Fixed underflow error on", underflow[1]
                else:
                    print "Problem!", underflow, probs[underflow[0]]
        

        allschemeprobs=map(lambda x:math.log(x, 2), allschemeprobs) 
        data_prob=sum(allschemeprobs)

        probs=e_norm_post(probs)  #normalize

        #M-step
        [t_table, rprobs]=m_frac_counts(words, stanzas, schemes, probs)
        
        #check convergence
        if ctr>0 and data_prob-old_data_prob<epsilon:
            break
        
        print "Iteration", ctr, "-- Log likelihood of data:", data_prob

        [t_table, rprobs]=m_norm_frac(t_table, len(words), rprobs)

    #error if it didn't converge
    if ctr==maxsteps-1 and data_prob-old_data_prob>=epsilon:
        print "Warning: EM did not converge"
    
    return [t_table, probs, data_prob]

def show_rhymes(probs, stanzas, schemes, outfile):
    """write rhyme schemes at convergence"""
    o=codecs.open(outfile, 'w', 'utf8')
    for stanza, stanzaprobs in zip(stanzas, probs):
        #scheme with highest probability
        bestscheme=schemes[len(stanza)][numpy.argmax(numpy.array(stanzaprobs))]
        o.write(' '.join(stanza)+'\n')
        o.write(' '.join(map(str, bestscheme))+'\n\n')
    o.close()

def init_uniform_r(schemes):
    """assign equal prob to every scheme"""
    rprobs={}
    numschemes=float(sum(map(len, schemes.values())))
    uni_prob=1/numschemes
    
    for i in schemes:
        for scheme in schemes[i]:
            rprobs[tuple(scheme)]=uni_prob

    return rprobs

def main(args):
    if len(args)!=3:
        print "Usage: findschemes.py gold-data init-type output-filename"
        print "where init-type may be u for uniform, o for orthographic"
        return
    
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  #to flush buffer
    
    #load stanzas and schemes
    INFILE=args[0]
    stanzas=load_stanzas(INFILE) 
    schemes=load_schemes('allschemes.pickle')
    print "Loaded files"
    
    #get list of words
    words=get_wordset(stanzas)

    #initialize p(r)
    rprobs=init_uniform_r(schemes)

    if args[1][0]=='u':  #uniform init
        t_table=init_uniform_ttable(words)
        
        print "Initialized,", len(words), "words"
        [final_t_table, final_probs, data_prob]=iterate(t_table, words, stanzas, schemes, rprobs, 100)
    
    elif args[1][0]=='o':  #init based on orthographic word sim
        t_table=init_basicortho_ttable(words)
        print "Initialized,", len(words), "words"
        [final_t_table, final_probs, data_prob]=iterate(t_table, words, stanzas, schemes, rprobs, 100)

    elif args[1][0]=='p':  #init based on rhyming definition
        t_table=init_perfect_ttable(words)
        print "Initialized,", len(words), "words"
        [final_t_table, final_probs, data_prob]=iterate(t_table, words, stanzas, schemes, rprobs, 100)
        
    #write rhyme schemes
    show_rhymes(final_probs, stanzas, schemes, args[2])
    print "Wrote result"

if __name__=='__main__':
    main(sys.argv[1:])
