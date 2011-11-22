#!/usr/bin/env python

"""Evaluate rhyme schemes against gold standard.
Also contains some utilities to parse data.
Jan 2011."""

import sys, pickle
import numpy, math, string
from collections import defaultdict
import codecs

def remove_punct(s):
    """remove non-letter chars"""
    return ''.join(filter(lambda c: c not in string.punctuation and c not in string.digits, s))

def get_wordset(poems):
    """get all words"""
    words=sorted(list(set(reduce(lambda x, y: x+y, poems))))
    return words

def get_rhymelists(poem, scheme):
    """transform poem into lists of rhymesets as given by rhyme scheme"""
    rhymelists=defaultdict(list)
    for poemword, schemeword in zip(poem, scheme):
        rhymelists[schemeword].append(poemword)
    return map(sorted, rhymelists.values())

patstr=list(string.lowercase+string.uppercase)+map(str, range(50))

def gen_pattern(seed, n):
    """generate scheme from seed"""
    seed=map(int, seed)
    if n%len(seed)>0:
        return "Error!"
    pattern=seed[:]
    increment=max(seed)
    while len(pattern)<n:
        pattern+=map(lambda seedunit: seedunit+increment, seed)
        increment=max(pattern)
    return map(str, pattern)

def parse(filename):
    """read rhyme schemes and poems/stanzas"""
    f=map(lambda x:x.split(), codecs.open(filename, encoding='utf-8').readlines())
    f=filter(lambda line: len(line)==0 or line[0] not in ['DATE', 'AUTHOR'], f)
    
    curscheme=[]
    curpoemscheme=[]
    stanzaschemes=[]
    poemschemes=[]
    curstanza=[]
    poems=[]
    stanzas=[]
    for i, line in enumerate(f):
        if len(line)==0:
            if curstanza!=[]:
                
                if curscheme==[]:
                    print "Error! No scheme read", curstanza
                
                stanzas.append(curstanza)
                stanzasize=len(curstanza)
                
                if curscheme[-1]=='*':
                    #generate pattern
                    genscheme=gen_pattern(curscheme[:-1], stanzasize)
                    if genscheme=='Error!':
                        print "Error! Seed doesn't match", i, curstanza, curscheme[:-1]
                    stanzaschemes.append(genscheme)
                elif len(curscheme)!=stanzasize:
                        print "Error! Stanza size and scheme size don't match", curstanza, curscheme
                else:
                    stanzaschemes.append(curscheme)

                if curpoemscheme[-1]=='*':
                    #generate pattern
                    genscheme=gen_pattern(curpoemscheme[:-1], stanzasize)
                    poemschemes.append(genscheme)
                elif len(curpoemscheme)!=stanzasize:
                        print "Error! Stanza size and scheme size don't match", curstanza, curpoemscheme
                else:
                    poemschemes.append(curpoemscheme)
                
                curstanza=[]
        
        elif line[0]=='RHYME':
            curscheme=map(lambda x: str(patstr.index(x)+1), line[1:-1])
            if line[-1]=='*':
                curscheme.append('*')
            else:
                curscheme.append(str(patstr.index(line[-1])+1))

            curpoemscheme=curscheme[:]  #in case RHYME-POEM isn't specified

        elif line[0]=='RHYME-POEM':
            curpoemscheme=map(lambda x: str(patstr.index(x)+1), line[1:])
        
        elif line[0]=='TITLE':
            if stanzas!=[]:
                poems.append(stanzas)
                stanzas=[]
        else:
            line=filter(lambda x:x!='', map(remove_punct, line))
            curstanza.append(line[-1].lower())

    if curstanza!=[]:
        stanzas.append(curstanza)
        stanzaschemes.append(curscheme)
        poemschemes.append(curpoemscheme)
    
    poems.append(stanzas)
    
    return [stanzaschemes, poemschemes, poems]

def dist_schemes(rhymeschemes, storefile, store=False):
    """all rhyme schemes of a given length, with frequencies"""
    dist=defaultdict(lambda : defaultdict(int))
    for r in rhymeschemes:
        dist[len(r)][' '.join(r)]+=1

    if not store:
        return dist

    newdist={}  #convert to normal dict with tuple keys in order to pickle
    for l in dist:
        newdist[l]=[]
        for r in dist[l]:
            newdist[l].append((r, dist[l][r]))
    allschemes=open(storefile, 'wb')
    pickle.dump(newdist, allschemes)
    allschemes.close()

    return dist

def save_gold_std(stanzaschemes, poemschemes, poems, filename):
    """write gold standard"""
    o=codecs.open(filename, 'w', 'utf-8')
    schemectr=0
    for pi, poem in enumerate(poems):
        for stanza in poem:
            stanzascheme=stanzaschemes[schemectr]
            poemscheme=poemschemes[schemectr]
            o.write('POEM'+str(pi)+' '+unicode(' '.join(stanza))+'\n')
            o.write(' '.join(stanzascheme)+'\n')
            o.write(' '.join(poemscheme)+'\n\n')
            schemectr+=1
    o.close()

def load_gold(filename):
    f=open(filename).readlines()
    stanzas=[]
    stanzaschemes=[]
    poemschemes=[]
    for i, line in enumerate(f):
        line=line.split()
        if i%4==0:
            stanzas.append(line[1:])
        elif i%4==1:
            if line==[]:
                print "Error in gold!", i, f[i-1], f[i-2]
            stanzaschemes.append(line)
        elif i%4==2:
            poemschemes.append(line)            
    return [stanzaschemes, poemschemes, stanzas]

def rhyming_entropy(stanzaschemes, stanzas):
    """compute entropy of rhyming pairs"""
    pairs=defaultdict(int)
    totalpairs=0.0
    for scheme, stanza in zip(stanzaschemes, stanzas):
        for i, (schemei, wordi) in enumerate(zip(scheme, stanza)):
            for (schemej, wordj) in zip(scheme[i+1:], stanza[i+1:]):
                totalpairs+=1
                if schemei==schemej:
                    if wordi<=wordj:
                        pairs[(wordi, wordj)]+=1
                    else:
                        pairs[(wordj, wordi)]+=1
    #normalize
    for pair in pairs:
        pairs[pair]=pairs[pair]/totalpairs
    #compute entropy
    return sum(map(lambda paircount: -paircount*math.log(paircount, 2), pairs.values()))

def scheme_entropy(stanzaschemes, stanzas):
    """compute entropy of rhyme schemes"""
    schemes=defaultdict(float)
    for scheme, stanza in zip(stanzaschemes, stanzas):
        schemes[tuple(scheme)]+=1.0
    #normalize
    total=len(stanzaschemes)
    for scheme in schemes:
        schemes[scheme]=schemes[scheme]/total
    #compute entropy
    return sum(map(lambda schemecount: -schemecount*math.log(schemecount, 2), schemes.values()))

def load_result(filename):
    f=open(filename).readlines()
    stanzas=[]
    schemes=[]
    for i, line in enumerate(f):
        line=line.split()
        if i%3==0:
            stanzas.append(line[1:])
        elif i%3==1:
            if line==[]:
                print "Error in result!", i, f[i-1], f[i-2]
            schemes.append(line)
    return [schemes, stanzas]

def stats(rhymeschemes):
    """show distribution of schemes of different lengths"""
    dist=defaultdict(lambda : defaultdict(int))
    for r in rhymeschemes:
        dist[len(r)][' '.join(r)]+=1
    for l in dist:
        print l, 
        distl=sorted(dist[l].items(), key=lambda x:x[1], reverse=True)
        for (r, v) in distl:
            print v,
        print

def compare(stanzas, gold_schemes, found_schemes):
    """get accuracy and precision/recall"""
    total=float(len(gold_schemes))
    correct=0.0
    for (g, f) in zip(gold_schemes, found_schemes):
        if g==f:
            correct+=1
    print "Accuracy", correct, total, 100*correct/total

    #for each word, let rhymeset[word] = set of words in rest of stanza rhyming with the word
    #precision = # correct words in rhymeset[word]/# words in proposed rhymeset[word]
    #recall = # correct words in rhymeset[word]/# words in reference words in rhymeset[word]
    #total precision and recall = avg over all words over all stanzas
    
    tot_p=0.0
    tot_r=0.0
    tot_words=0.0
    
    for (s, g, f) in zip(stanzas, gold_schemes, found_schemes):
        stanzasize=len(s)
        for wi, word in enumerate(s):
            grhymeset_word = set(map(lambda x:x[0], filter(lambda x:x[1]==g[wi], zip(range(wi+1, stanzasize), g[wi+1:]))))
            frhymeset_word = set(map(lambda x:x[0], filter(lambda x:x[1]==f[wi], zip(range(wi+1, stanzasize), f[wi+1:]))))

            if len(grhymeset_word)==0:
                continue

            tot_words+=1

            if len(frhymeset_word)==0:
                continue
            
            #find intersection
            correct=float(len(grhymeset_word.intersection(frhymeset_word)))
            precision=correct/len(frhymeset_word)
            recall=correct/len(grhymeset_word)
            tot_p+=precision
            tot_r+=recall

    precision=tot_p/tot_words
    recall=tot_r/tot_words
    print "Precision", precision
    print "Recall", recall
    print "F-score", 2*precision*recall/(precision+recall)
    
def naive(gold_schemes):
    """find naive baseline (most common scheme of a given length)?"""
    dist = pickle.load(open('allschemes.pickle'))
    best_schemes={}
    for i in dist:
        if dist[i]==[]:
            continue
        best_schemes[i]=(max(dist[i], key=lambda x:x[1])[0]).split()

    naive_schemes=[]
    for g in gold_schemes:
        naive_schemes.append(best_schemes[len(g)])
    return naive_schemes

def lessnaive(gold_schemes):
    """find 'less naive' baseline (most common scheme of a given length in subcorpus)"""
    best_schemes=defaultdict(lambda : defaultdict(int))
    for g in gold_schemes:
        best_schemes[len(g)][tuple(g)]+=1

    m=sum(map(len, best_schemes.values()))
    
    for i in best_schemes:
        best_schemes[i]=list(max(best_schemes[i].items(), key=lambda x:x[1])[0])

    naive_schemes=[]
    for g in gold_schemes:
        naive_schemes.append(best_schemes[len(g)])
    return naive_schemes

def main(args):
    if len(args)<1 or len(args)>2:
        print "Usage: evaluate.py gold-file [hypothesis-output-filename]"
        return 
    
    GOLD=args[0]    
    [gstanzaschemes, gpoemschemes, gstanzas]=load_gold(GOLD)

    words=get_wordset(gstanzas)
    n=len(words)
    
    #for stanzas 
    print 'Num of stanzas: ', len(gstanzas)
    print 'Num of lines: ', sum(map(len, gstanzas))
    print 'Num of end word types: ', len(words)
    print
    
    naive_schemes=naive(gstanzaschemes)
    print "Naive baseline:"
    compare(gstanzas, gstanzaschemes, naive_schemes)
    print

    lessnaive_schemes=lessnaive(gstanzaschemes)
    print "Less naive baseline:"
    compare(gstanzas, gstanzaschemes, lessnaive_schemes)
    print

    if len(args)>1:
        HYP=args[1]
        [hstanzaschemes, hstanzas]=load_result(HYP)
        print HYP,":"
        compare(gstanzas, gstanzaschemes, hstanzaschemes)
        print

if __name__=='__main__':
    main(sys.argv[1:])
