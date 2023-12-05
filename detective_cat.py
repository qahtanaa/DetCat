import sys
from functools import lru_cache
import pandas as pd
import numpy as np
import csv
import scipy
from scipy.cluster.hierarchy import fclusterdata
import itertools
from scipy import stats
from typing_extensions import Counter
import ntpath
from functools import reduce
from strsimpy.string_distance import StringDistance
from strsimpy.weighted_levenshtein import WeightedLevenshtein


# UL=["A",'B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
UL = list(map(chr, range(ord('A'), ord('Z')+1)))
LL = list(map(chr, range(ord('a'), ord('z')+1)))
# DD = ["0",'1','2','3','4','5','6','7','8','9']
DD = list(map(chr, range(ord('0'), ord('9')+1)))
S = list()
WS = [" ", "\t", "\n"]

def initialize_S(word_list):
  global S
  for word in word_list:
    for c in word:
      if (not(c in UL) and not(c in LL) and not(c in DD) and not (c in WS) and not(c in S)):
        S.append(c)

def abstractwords(word): #this function builds and returns a list that contains a list of each value's literal character and their class abstraction. [char,abstraction]
  wordlist=[]
  #word=word.strip('"') done for results purposes only for rayyan file (also done for rayyan results)
  for char in word:
    if char in UL:
      wordlist.append([char,"Ů"])
    elif char in LL:
      wordlist.append([char,"Ł"])
    elif char in DD:
      wordlist.append([char,"Đ"])
    elif char in WS:
      wordlist.append([char])
    else:
      wordlist.append([char,"Š"])
  return(wordlist)


coverage={}
g=[]

def read_table(tab_name):
    t_name = ntpath.basename(tab_name)
    try:
        df = pd.read_csv(filepath_or_buffer=tab_name, dtype=object, delimiter=',', low_memory=False,
                         quoting=csv.QUOTE_ALL, doublequote=True)
    except ValueError:
        try:
            df = pd.read_csv(filepath_or_buffer=tab_name, dtype=object, delimiter=',', low_memory=False,
                             quoting=csv.QUOTE_ALL, doublequote=True, encoding="ISO-8859-1")
        except:
            print("Error reading csv file .. file encoding is not recognizable")
            return None
    return df

def getpatterns(d,word,limlength,maxlen): #this function returns generator of patterns
    #word=word.strip('"') #added for journal name. unreadable output so removing the " sign for both.
    if len(word)<limlength:
      for i in word:
        coverage.update({word:''})
      g.clear()
      combinations = itertools.product(*d)
      flat_list = [item for sublist in d for item in sublist]
      flat_list= (''.join(flat_list))
      for c in combinations:
          c=''.join(c)
          g.append(c)
      coverage.update({word:g[:]})
      return(g)
    elif len(word)>limlength and len(word)<maxlen:
      a=word.split() #splitting by words.
      longword=[]
      for i in a:
        j = i
        for char in i:
          if char in UL:     # Uppercase letters
            j = j.replace(char,"Ů")
          elif char in LL:   # Lowercase letters
            j = j.replace(char,"Ł")
          elif char in DD:   # digits
            j = j.replace(char,"Đ")
          else:       # Everything else
            j = j.replace(char,"Š")
          #elif char in WS:  chose not to abstract white spaces because it's just one value so creating redundant patterns.
            #j.replace(char,"⊔")
        longword.append([i,j])
      g.clear()
      combinations = itertools.product(*longword)
      for c in combinations:
        c=(' '.join(c))
        g.append(c)
      coverage.update({word:g[:]})
      return(g)

    else:
      verylongword=[]
      g.clear()
      j=word
      for char in word:
        if char in UL:     # Uppercase letters
          j=j.replace(char,"Ů")
        elif char in LL:   # Lowercase letters
          j=j.replace(char,"Ł")
        elif char in DD:   # digits
          j=j.replace(char,"Đ")
        else:       # Everything else
          j=j.replace(char,"Š")
      #print(word,j)
      g.append(word)
      g.append(j)
      coverage.update({word:g[:]})
      return(g)



b={}
def getdictionary(x,count_of_word): #this function checks if a pattern exists in the dictionary, and if it does not, creates a new entry. if it does exist, adds to counter.
  for i in x:
    counter=b.get(i,count_of_word) # in the beginning, b is an empty dict. counter checks if i, a pattern in list/dict x is in b, if it is, we assign the current freq to the counter, otherwise, we assign 1 to a new key.
    if i in b:
      counter=counter+count_of_word
    b.update({i:counter})
  return(b)

def getscore(b, genlevelcostmultiplier, minpatmultiplier, acc_threshold, 
            lencomparison,pvalue,smallnfpvalue, smallnfnvalue):
  sorted_dict = dict(sorted(b.items(), key=lambda x:x[1],reverse=True))
  sorted_dict2=sorted_dict.items()
  #genlevelcostmultiplier=1.2
  for i in sorted_dict.items():
    genlevel=0
    for char in i[0]:
      if char in ('ŮŁĐŠ⊔'):
        genlevel=genlevel+1

    score=((i[1])*(1/((genlevel*genlevelcostmultiplier)+1)))
    sorted_dict.update({i[0]:score})
  sorted_dict = dict(sorted(sorted_dict.items(), key=lambda x:x[1],reverse=True))
  return(loopthrough(sorted_dict,minpatmultiplier,acc_threshold,lencomparison,pvalue,smallnfpvalue, smallnfnvalue))


def loopthrough(sorted_dict,minpatmultiplier,acc_threshold,lencomparison,pvalue,smallnfpvalue, smallnfnvalue):
  def get_rows(i, patternlist):
    coverage_values=list(coverage.values())
    subtractthis=[]
    saveforlater=[]
    for row in coverage_values:
      if patternlist[i] in row:
        subtractthis.append([x for x in row if x != patternlist[i]])
        saveforlater.append(row[0])
    pattern_values.update({patternlist[i]:saveforlater})
    return subtractthis

  def get_outliers(patternlist,minpatmultiplier,acc_threshold,lencomparison,pvalue,smallnfpvalue, smallnfnvalue):
    freq_cov={}
    outlier_patterns=[]
    dominant_patterns=[]
    valid_values=[]
    outlier_values=[]
    final_patterns={}
    final_values={}
    accumulate=0
    #acc_threshold=0.95
    for i in patternlist:
      freq_cov.update({i:b.get(i)})
    total=sum(freq_cov.values())
    freq_cov_2=list(freq_cov.values())
    mincov=(np.mean(freq_cov_2)/total)*minpatmultiplier
    for value in freq_cov.items():
      coveragevalue=(value[1]/total)
      if coveragevalue > mincov and accumulate < acc_threshold:
        dominant_patterns.append(value[0])
        accumulate=accumulate+coveragevalue
      elif (value[1] > 50):
        dominant_patterns.append(value[0])
        accumulate = accumulate+coveragevalue
      else:
        outlier_patterns.append(value[0])
      # print(coveragevalue, accumulate)
    for dom_pat in dominant_patterns:
      valid_values.append(pattern_values.get(dom_pat))
    for out_pat in outlier_patterns:
      outlier_values.append(pattern_values.get(out_pat))
    valid_values = [item for sublist in valid_values for item in sublist]
    outlier_values= [item for sublist in outlier_values for item in sublist]
    final_patterns.update({'dominant': dominant_patterns})
    final_patterns.update({'outlier patterns': outlier_patterns})
    final_values.update({'valid':valid_values})
    final_values.update({'outliers':outlier_values})

    return(get_distances(freq_cov,patternlist,outlier_patterns,dominant_patterns, final_values,final_patterns,lencomparison,pvalue,smallnfpvalue, smallnfnvalue))


  def get_distances(freq_cov,patternlist,outlier_patterns,dominant_patterns,final_values, final_patterns,lencomparison,pvalue,smallnfpvalue, smallnfnvalue): #calculates distances using class-weighted measure between pairwise dominant-dominant... and dominant-outlier. performs statistical z test to determine if significant results.    
    domdom=itertools.product(dominant_patterns,dominant_patterns)
    domout=itertools.product(dominant_patterns,outlier_patterns)
    domoutscores={}
    domdomscores={}
    fpposition=[]
    fppositiondomdom=[]
    false_positives=[]
    false_negatives=[]
    for i in domout:
      if (len(i[0])<lencomparison) and (len(i[1])<lencomparison): #similarity measure length limitation
        try:
          domoutscores.update({i:AZWeightedLevenshtein.distance(i[0],i[1])})
        except:
          if len(i[0]) != len(i[0].encode()):
            print(i, "cannot handle this character")
    for i in domdom:
      if (len(i[0])<lencomparison) and (len(i[1])<lencomparison): #limiting pairwise distance calculation for efficiency
        try:
          domdomscores.update({i:AZWeightedLevenshtein.distance(i[0],i[1])})
        except:
          if len(i[0]) != len(i[0].encode()):
            print(i, "cannot handle this character")
    domoutscoreslist=list(domoutscores.values())
    domoutarr = np.array(domoutscoreslist)
    domoutz=stats.zscore(domoutarr)
    if len(domoutz)>15:
      for index,i in enumerate(domoutz):
        if i <0 and scipy.stats.norm.sf(abs(i)) <pvalue:
          fpposition.append(index)
      for i in fpposition:
        if list(domoutscores.keys())[i][1] not in false_positives:
          false_positives.append(list(domoutscores.keys())[i][1])
    else:
       for item in domoutscores.items():
         if (item[1])<smallnfpvalue and item[0][1] not in false_positives:
           false_positives.append(item[0][1])
    domdomscoreslist=list(domdomscores.values())
    domdomarr = np.array(domdomscoreslist)
    domdomz=stats.zscore(domdomarr)
    if len(domdomz)>15:
      for index,i in enumerate(domdomz):
       if i >0 and scipy.stats.norm.sf(abs(i)) <pvalue:
        fppositiondomdom.append(index)
      for i in fppositiondomdom:
        if (freq_cov[list(domdomscores.keys())[i][0]]>freq_cov[list(domdomscores.keys())[i][1]]):
          if list(domdomscores.keys())[i][1] not in false_negatives:
           false_negatives.append(list(domdomscores.keys())[i][1])
        else:
          if list(domdomscores.keys())[i][0] not in false_negatives:
            false_negatives.append(list(domdomscores.keys())[i][0])
    else:
      for item in domdomscores.items():   #this part is added in case the number of pairwise comparisons is smaller than 15, making statistical tests less valid.
        if(item[1]>smallnfnvalue): #largest cost of substitution of one inter-class character.
          if (freq_cov.get(item[0][0])>freq_cov.get(item[0][1])):
            if item[0][1] not in false_negatives:
              false_negatives.append(item[0][1])
          else:
            if item[0][0] not in false_negatives:
              false_negatives.append(item[0][0])
    fp_values=[]
    fn_values=[]
    for fp_pat in false_positives:
        fp_values.append(pattern_values.get(fp_pat))
    for fn_pat in false_negatives:
        fn_values.append(pattern_values.get(fn_pat))
    fp_values= [item for sublist in fp_values for item in sublist]
    fn_values= [item for sublist in fn_values for item in sublist]
    final_values.update({"FP":fp_values})
    final_values.update({"FN":fn_values})
    final_patterns.update({"domdomscores: ": domdomscores})
    final_patterns.update({"domoutscores: ": domoutscores})
    final_patterns.update({"FP_pat": false_positives})
    final_patterns.update({"FN_pat": false_negatives})
    df = pd.DataFrame.from_dict(final_values,orient='index').transpose()
    # print('values:',final_values, '\npatterns',final_patterns)
    
    df.to_csv('results_col.csv', index=False)
    # print('*'*10)
    # print(freq_cov)
    
    output_patterns = dict()
    output_patterns.clear()
    output_patterns['dominant'] = dict()
    output_patterns['outlier patterns'] = dict()
    for k in freq_cov.keys():
      if (k in final_patterns['dominant']):
        output_patterns['dominant'][k] = freq_cov[k]
      if (k in final_patterns['outlier patterns']):
        output_patterns['outlier patterns'][k] = freq_cov[k]
    # print('=='*24)
    # print('Domenant Patterns', final_patterns['dominant'])
    # print('=='*24)
    # print('Outlier Patterns: ', final_patterns['outlier patterns'])
    # print('=='*24)
    # df.to_csv('results_col.csv', index=False)
    # print('*'*10)
    # print(freq_cov)
    # print('*'*10)
    return output_patterns

  coverage_values=list(coverage.values())
  patternlist=list(sorted_dict.keys())
  plistset=set(patternlist)
  pattern_values={}
  for i in range (len(sorted_dict)):
    if i <len(patternlist):
      remainingrow=get_rows(i,patternlist)
      remainingrowflat = [item for sublist in remainingrow for item in sublist]
      for m in remainingrowflat:
        if m in plistset:
          patternlist.remove(m)
          plistset.remove(m)
    else:
      break
  # print('+'*10)
  # print('patterns = ', patternlist)
  # print('+'*10)
  return(get_outliers(patternlist,minpatmultiplier,acc_threshold,lencomparison,pvalue,smallnfpvalue, smallnfnvalue))


def default_insertion_cost(char):
    return 1.0
def default_deletion_cost(char):
    return 1.0
def default_substitution_cost(char_a, char_b):
    return 1.0

class AZWeightedLevenshtein(StringDistance):
    def __init__(self,
                 substitution_cost_fn=default_substitution_cost,
                 insertion_cost_fn=default_insertion_cost,
                 deletion_cost_fn=default_deletion_cost,
                 ):
        self.substitution_cost_fn = substitution_cost_fn
        self.insertion_cost_fn = insertion_cost_fn
        self.deletion_cost_fn = deletion_cost_fn

    def distance(self, s0, s1):
        if s0 is None:
            raise TypeError("Argument s0 is NoneType.")
        if s1 is None:
            raise TypeError("Argument s1 is NoneType.")
        if s0 == s1:
            return 0.0
        if len(s0) == 0:
            return reduce(lambda cost, char: cost + self.insertion_cost_fn(char,s0), s1, 0)
        if len(s1) == 0:
            return reduce(lambda cost, char: cost + self.deletion_cost_fn(char,s0), s0, 0)

        v0, v1 = [0.0] * (len(s1) + 1), [0.0] * (len(s1) + 1)
        v0[0] = 0
        for i in range(1, len(v0)):
            v0[i] = v0[i - 1] + self.insertion_cost_fn(s1[i - 1],s0)

        for i in range(len(s0)):
            s0i = s0[i]
            deletion_cost = self.deletion_cost_fn(s0i,s0)
            v1[0] = v0[0] + deletion_cost

            for j in range(len(s1)):
                s1j = s1[j]
                cost = 0
                if s0i != s1j:
                    cost = self.substitution_cost_fn(s0i, s1j)
                insertion_cost = self.insertion_cost_fn(s1j,s0)
                v1[j + 1] = min(v1[j] + insertion_cost, v0[j + 1] + deletion_cost, v0[j] + cost)
            v0, v1 = v1, v0

        return v0[len(s1)]


def insertion_cost(char,word):
    if (word[len(word)-1] in LL or word[len(word)-1]=='Ł')  and (char in LL or char=="Ł"):
      return 1.5
    elif (word[len(word)-1] in LL or word[len(word)-1]=='Ł') and (char not in LL or char !="Ł"):
      return 6
    if (word[len(word)-1] in UL or word[len(word)-1]=='Ů') and (char in UL or char=="Ů"):
      return 1.5
    elif (word[len(word)-1] in UL or word[len(word)-1]=='Ů') and (char not in UL or char !="Ů"):
      return 6
    if (word[len(word)-1] in DD or word[len(word)-1]=='Đ') and (char in DD or char=="Đ"):
      return 1.5
    elif (word[len(word)-1] in DD or word[len(word)-1]=='Đ') and (char not in DD or char !="Đ"):
      return 6
    if (word[len(word)-1] in S or word[len(word)-1]=='Š') and (char in S or char=="Š"):
      return 1.5
    elif (word[len(word)-1] in S or word[len(word)-1]=='Š') and (char not in S or char !="Š"):
      return 6
    if (word[len(word)-1] in WS or word[len(word)-1]=='⊔') and (char in WS or char=="⊔"):
      return 1.5
    elif (word[len(word)-1] in WS or word[len(word)-1]=='⊔') and (char not in WS or char !="⊔"):
      return 6

def deletion_cost(char,word):
    if (word[len(word)-2] in LL or word[len(word)-2]=='Ł')  and (char in LL or char=="Ł"):
      return 1.5
    elif (word[len(word)-2] in LL or word[len(word)-2]=='Ł') and (char not in LL or char !="Ł"):
      return 6
    if (word[len(word)-2] in UL or word[len(word)-2]=='Ů') and (char in UL or char=="Ů"):
      return 1.5
    elif (word[len(word)-2] in UL or word[len(word)-2]=='Ů') and (char not in UL or char !="Ů"):
      return 6
    if (word[len(word)-2] in DD or word[len(word)-2]=='Đ') and (char in DD or char=="Đ"):
      return 1.5
    elif (word[len(word)-2] in DD or word[len(word)-2]=='Đ') and (char not in DD or char !="Đ"):
      return 6
    if (word[len(word)-2] in S or word[len(word)-2]=='Š') and (char in S or char=="Š"):
      return 1.5
    elif (word[len(word)-2] in S or word[len(word)-2]=='Š') and (char not in S or char !="Š"):
      return 6
    if (word[len(word)-2] in WS or word[len(word)-2]=='⊔') and (char in WS or char=="⊔"):
      return 1.5
    elif (word[len(word)-2] in WS or word[len(word)-2]=='⊔') and (char not in WS or char !="⊔"):
      return 6


def substitution_cost(char_a, char_b):
    if (char_a == "Ł" and char_b in LL) or (char_a in LL and char_b =="Ł"):
        return 0.5
    elif (char_a == "Ł" and char_b in UL) or (char_a in UL and char_b =="Ł"):
        return 2.5
    elif (char_a == "Ł" and char_b =="Ů") or (char_a =="Ů" and char_b =="Ł"):
        return 2
    elif (char_a == "Ł" and char_b in DD) or (char_a in DD and char_b =="Ł"):
        return 7.5
    elif (char_a == "Ł" and char_b =="Đ") or (char_a =="Đ" and char_b =="Ł"):
        return 7
    elif (char_a == "Ł" and char_b in S) or (char_a in S and char_b =="Ł"):
        return 7.5
    elif (char_a == "Ł" and char_b =="Š") or (char_a =="Š" and char_b =="Ł"):
        return 7
    elif (char_a == "Ł" and char_b in WS) or (char_a in WS and char_b =="Ł"):
        return 7.5
    elif (char_a == "Ł" and char_b =="⊔") or (char_a =="⊔" and char_b =="Ł"):
        return 7

    if (char_a == "Ů" and char_b in UL) or (char_a in UL and char_b =="Ů"):
        return 0.5
    elif (char_a == "Ů" and char_b in LL) or (char_a in LL and char_b =="Ů"):
        return 2.5
    elif (char_a == "Ů" and char_b in DD) or (char_a in DD and char_b =="Ů"):
        return 7.5
    elif (char_a == "Ů" and char_b =="Đ") or (char_a =="Đ" and char_b =="Ů"):
        return 7
    elif (char_a == "Ů" and char_b in S) or (char_a in S and char_b =="Ů"):
        return 7.5
    elif (char_a == "Ů" and char_b =="Š") or (char_a =="Š" and char_b =="Ů"):
        return 7
    elif (char_a == "Ů" and char_b in WS) or (char_a in WS and char_b =="Ů"):
        return 7.5
    elif (char_a == "Ů" and char_b =="⊔") or (char_a =="⊔" and char_b =="Ů"):
        return 7

    if (char_a == "Đ" and char_b in DD) or (char_a in DD and char_b =="Đ"):
        return 0.5
    elif (char_a == "Đ" and char_b in LL) or (char_a in LL and char_b =="Đ"):
        return 7.5
    elif (char_a == "Đ" and char_b in UL) or (char_a in UL and char_b =="Đ"):
        return 7.5
    elif (char_a == "Đ" and char_b in S) or (char_a in S and char_b =="Đ"):
        return 6.5
    elif (char_a == "Đ" and char_b =="Š") or (char_a =="Š" and char_b =="Đ"):
        return 6
    elif (char_a == "Đ" and char_b in WS) or (char_a in WS and char_b =="Đ"):
        return 6.5
    elif (char_a == "Đ" and char_b =="⊔") or (char_a =="⊔" and char_b =="Đ"):
        return 6

    if (char_a == "Š" and char_b in S) or (char_a in S and char_b =="Š"):
        return 0.5
    elif (char_a == "Š" and char_b in LL) or (char_a in LL and char_b =="Š"):
        return 7.5
    elif (char_a == "Š" and char_b in UL) or (char_a in UL and char_b =="Š"):
        return 7.5
    elif (char_a == "Š" and char_b in DD) or (char_a in DD and char_b =="Š"):
        return 6.5
    elif (char_a == "Š" and char_b in WS) or (char_a in WS and char_b =="Š"):
        return 6.5
    elif (char_a == "Š" and char_b =="⊔") or (char_a =="⊔" and char_b =="Š"):
        return 6

    if (char_a == "⊔" and char_b in WS) or (char_a in WS and char_b =="⊔"):
        return 0.5
    elif (char_a == "⊔" and char_b in LL) or (char_a in LL and char_b =="⊔"):
        return 7.5
    elif (char_a == "⊔" and char_b in UL) or (char_a in UL and char_b =="⊔"):
        return 7.5
    elif (char_a == "⊔" and char_b in DD) or (char_a in DD and char_b =="⊔"):
        return 6.5
    elif (char_a == "⊔" and char_b in S) or (char_a in S and char_b =="⊔"):
        return 6.5

#Literals:
    if (char_a in UL and char_b in LL) or (char_a in LL and char_b in UL):
        return 3
    elif (char_a in UL and char_b in DD) or (char_a in DD and char_b in UL):
        return 8
    elif (char_a in UL and char_b in S) or (char_a in S and char_b in UL):
        return 8
    elif (char_a in UL and char_b in WS) or (char_a in WS and char_b in UL):
        return 8

    if (char_a in LL and char_b in DD) or (char_a in DD and char_b in LL):
        return 8
    elif (char_a in LL and char_b in S) or (char_a in S and char_b in LL):
        return 8
    elif (char_a in LL and char_b in WS) or (char_a in WS and char_b in LL):
        return 8

    if (char_a in DD and char_b in S) or (char_a in S and char_b in DD):
        return 7
    elif (char_a in DD and char_b in WS) or (char_a in WS and char_b in DD):
        return 7

    if (char_a in S and char_b in WS) or (char_a in WS and char_b in S):
        return 7

    if (char_a in UL and char_b in UL) or (char_a in LL and char_b in LL) or (char_a in DD and char_b in DD) or (char_a in WS and char_b in WS) or (char_a in S and char_b in S):
        return 0.5
    return 1.0

AZWeightedLevenshtein = AZWeightedLevenshtein(
    substitution_cost_fn=substitution_cost,
    insertion_cost_fn=insertion_cost,
    deletion_cost_fn=deletion_cost)



def get_att_details(att):
    avg_len = 0
    max_len = 0
    num_int = 0
    num_float = 0
    num_text = 0
    thresh = 0.99
    details = dict()
    dtype_dict = dict()
    dtype_dict["typeInt"] = dict()
    dtype_dict["typeFloat"] = dict()
    dtype_dict["typeText"] = dict()
    dtype = ""
    no_non_empty = 0
    att_no_null_values = att.dropna()
    idxA = att_no_null_values.index.tolist()
    no_non_empty = len(att_no_null_values)
    min_len = len(att_no_null_values[idxA[0]])
    for ix in idxA:
        # find min, max, average length of the values
        v = att_no_null_values[ix]
        lv = len(str(v))
        avg_len += lv
        max_len = max(max_len, lv)
        min_len = min(min_len, lv)

        # check the types of the values
        try:
            val = int(v)
            num_int += 1
            if lv in dtype_dict["typeInt"]:
                dtype_dict["typeInt"][lv] += 1
            else:
                dtype_dict["typeInt"][lv] = 1
        except ValueError:
            try:
                val = float(v)
                num_float += 1
                if lv in dtype_dict["typeFloat"]:
                    dtype_dict["typeFloat"][lv] += 1
                else:
                    dtype_dict["typeFloat"][lv] = 1
            except ValueError:
                num_text += 1
                if lv in dtype_dict["typeText"]:
                    dtype_dict["typeText"][lv] += 1
                else:
                    dtype_dict["typeText"][lv] = 1
    if no_non_empty > 0:
        avg_len /= no_non_empty
    details["max_len"] = max_len
    details["avg_len"] = avg_len
    details["min_len"] = min_len
    # print(num_int, no_non_empty, thresh, num_int / no_non_empty)
    if num_int / no_non_empty > thresh:
        lens = sorted(dtype_dict["typeInt"].items(), key=lambda kv: kv[1], reverse=True)
        aa, bb = lens[0]
        if (bb / num_int > thresh) and (min_len >= 2):
            dtype = "code"
        else:
            dtype = "integer"        
        # print(lens) 
    if num_float / no_non_empty > thresh:
        dtype = "float"
    if dtype == "":
        lens = sorted(dtype_dict["typeText"].items(), key=lambda kv: kv[1], reverse=True)
        if len(lens) > 0:
            aa, bb = lens[0]
            if bb / num_text > thresh:
                dtype = "code"
            else:
                dtype = "text"
            details["lens"] = lens
        else:
            details["lens"] = []
    details["dtype"] = dtype
    
    return details


def get_df_details(df):
    df_details = dict()
    df_details.clear()
    for att_id in range(len(df.columns)):
        att = df[df.columns[att_id]]
        details = get_att_details(att)
        df_details[att_id] = details
        df_details[att_id]["att_name"] = df.columns[att_id]
    return df_details




# def main(df, limlength,maxlen,genlevelcostmultiplier, minpatmultiplier,acc_threshold,lencomparison,pvalue,smallnfpvalue, smallnfnvalue):
  # df = pd.read_csv(r'data/beers/beersDirty.csv',sep=',')
  # df = df.astype(str)
  
def find_cand_atts(df_dets):
  cand_atts = []
  for cc in range(len(cand_atts)):
        cand_atts.remove(cand_atts[0])
  for att_id in df_dets.keys():
    if (df_dets[att_id]["dtype"] == "integer") or (df_dets[att_id]["dtype"] == "float"):
      continue
    cand_atts.append(att_id)
  return cand_atts


def run_dc(param_config):
  global S 
  tab_name = param_config["tab_name"]
  df = read_table(tab_name)
  acc_threshold = param_config["min_acceptable_coverage"]
  limlength = param_config["MaxLen"]
  maxlen = 40
  genlevelcostmultiplier = 1.2
  minpatmultiplier = 0.05
  
  lencomparison = 20
  pvalue = 0.05
  smallnfpvalue = 1.5
  smallnfnvalue = 8
  df_dets = get_df_details(df)
  final_results = dict()
  # for k in df_dets.keys():
  #   print(k, df_dets[k]["dtype"])
  # return 
  cand_list = find_cand_atts(df_dets)
  for att in cand_list:
      cur_att = df.columns[att]
      # if cur_att == 'SSN':
      #   continue
      print('='*64)
      print("Working on :", cur_att)
      print('='*64)
      testlist = pd.Series.to_dict((df[cur_att].value_counts()))
      initialize_S(testlist)
      print(S)

      b.clear()
      g.clear()
      coverage.clear()

      for word in testlist.items():
        wordlist=abstractwords(word[0])
        patterns=getpatterns(wordlist,word[0],limlength,maxlen)
        getdict=getdictionary(patterns,word[1])
      # print("# of total patterns:",len(b))
      patterns_and_outliers = getscore(b, genlevelcostmultiplier, minpatmultiplier, 
                acc_threshold, lencomparison, pvalue, smallnfpvalue, smallnfnvalue)
      # print('='*64)
      # print(final)
      final_results[cur_att] = patterns_and_outliers
      final_results[cur_att]["dType"] = df_dets[att]["dtype"]
      # print('='*64)
      # print(cur_att)
      # print('='*64)
      # print(patterns_and_outliers)
      S.clear()
      # break
  print('Done .. ')
  return final_results
  
  


def run_dc_local(tab_name):
  try:
    # T = pd.read_csv(table_name, dtype=str, keep_default_na=False)
    df = read_table(tab_name)
  except OSError as e:
    print("Error reading csv!")
    sys.exit(1)
  acc_threshold = 0.95
  limlength = 15
  maxlen = 40
  genlevelcostmultiplier = 1.2
  minpatmultiplier = 0.7  
  lencomparison = 20
  pvalue = 0.05
  smallnfpvalue = 1.5
  smallnfnvalue = 8
  df_dets = get_df_details(df)
  final_results = dict()

  for k in df_dets.keys():
    print(k, df.columns[k], df_dets[k]["dtype"])    
    print('*'*20)
    print (df_dets[k])
    print('*'*20)
  return 
  cand_list = find_cand_atts(df_dets)
  for att in cand_list:
      cur_att = df.columns[att]
      print('='*64)
      print("Working on :", cur_att)
      print('='*64)
      testlist = pd.Series.to_dict((df[cur_att].value_counts()))
      initialize_S(testlist)
      print(S)

      b.clear()
      g.clear()
      coverage.clear()

      for word in testlist.items():
        wordlist=abstractwords(word[0])
        patterns=getpatterns(wordlist,word[0],limlength,maxlen)
        getdict=getdictionary(patterns,word[1])
      # print("# of total patterns:",len(b))
      patterns_and_outliers = getscore(b, genlevelcostmultiplier, minpatmultiplier, 
                acc_threshold, lencomparison, pvalue, smallnfpvalue, smallnfnvalue)
      # print('='*64)
      # print(final)
      final_results[cur_att] = patterns_and_outliers
      S.clear()
      # break
  print(final_results)
  


def main():
    # check arguments
    # for arg in sys.argv[1:]:
    if(len(sys.argv) != 2):
        print("Wrong number of arguments .. entered (",len(sys.argv),")")
        # print(sys.argv, file=sys.stderr)
        print("Usage (",sys.argv[0],"): <data file name>")
        sys.exit(1)

    table_name = sys.argv[1];
    
    run_dc_local(table_name)
    





if __name__ == "__main__":
    main()