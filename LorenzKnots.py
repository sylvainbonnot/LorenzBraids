from collections import OrderedDict
import numpy as np


def start(pair):
    return pair[0]

def end(pair):
    return pair[1]

def crossing_symbol(pair, symbol):
    crossings = []
    if symbol=='LL':
        for item in list_RL:
            if end(item)<end(pair):
                crossings.append(1)
    if symbol=='LR':
        for item in list_RL:
            crossings.append(1)
        for item in list_RR:
            if end(item)<end(pair):
                crossings.append(1)
    if symbol=='RL':
        for item in list_LR[::-1]:
            crossings.append(0)
        for item in list_LL[::-1]:
            if end(item)>end(pair):
                crossings.append(0)
                
    if symbol=='RR':
        for item in list_LR[::-1]:
            if end(item)>end(pair):
                crossings.append(0)
                
    return crossings

def word_to_cyclic_permutation(word):
    res = []
    for i in range(len(word)):
        res.append(word[i:]+word[:i])
    return sorted(res)

def permuted_list(word):
    lex_ordered_list = word_to_cyclic_permutation(word)
    res = []
    for item in lex_ordered_list:
        res.append(cyclic_permute(item))
    return res

def cyclic_permute(word):
    return word[1:]+word[:1]

def reverse_dict(dic):
    return dict(((value,key) for (key,value) in dic.items()))

def pairs_crossing(word):
    list_word = word_to_cyclic_permutation(word)
    list_word2 = permuted_list(word)
    list_LL =[]
    list_LR = []
    list_RL = []
    list_RR = []
    #res =[]
    for a,b in zip(list_word, list_word2):
        if a[0]=='L' and b[0]=='L':
            list_LL.append((a,b))
        elif a[0]=='L' and b[0]=='R':
            list_LR.append((a,b))
        elif a[0]=='R' and b[0]=='L':
            list_RL.append((a,b))
        elif a[0]=='R' and b[0]=='R':
            list_RR.append((a,b))
    return list_LL, list_LR, list_RL, list_RR

def all_symbols_pairs(word):
    double_word = word*2
    pairs = [double_word[i:i+2] for i in range(len(word))]
    return pairs

def odd_order(pair):
    a, b =pair
    
    if abs(a)%2==1:
        return a,b
    else:
        return b,a

class LorenzKnot():
    def __init__(self, word):
        self.word = word
        self.first_word = word_to_cyclic_permutation(self.word)[0]
        self.word_length = len(word)
        self.list_LL, self.list_LR, self.list_RL, self.list_RR=pairs_crossing(self.first_word)
        self.ll_num, self.lr_num, self.rl_num, self.rr_num = self.get_list_numerical_pairs()
              
        
        
    def list_lex_permutations(self):
        return word_to_cyclic_permutation(self.word)
    
    def list_cyclic_permutations(self):
        res = []
        current = self.first_word
        for i in range(self.word_length):
            res.append(current)
            current = cyclic_permute(current)
        return res
    
    def dict_lex_permutations(self):
        list_words = word_to_cyclic_permutation(self.word)
        return dict(zip(list_words, range(len(list_words))))
    
    def rev_dict_lex_permutations(self):
        return reverse_dict(self.dict_lex_permutations())
    
    def induced_permutation(self):
        dd = self.dict_lex_permutations()
        return [dd[cyclic_permute(item)] for item  in self.list_lex_permutations()]
    
    def integer_pairs(self):
        return zip(range(self.word_length), self.induced_permutation())
    
    def path_knot(self):
        list_cyclic = self.list_cyclic_permutations()
        dd = self.dict_lex_permutations()

        return [dd[item] for item in list_cyclic]
     
    def list_strands(self):
        path = self.path_knot()
        path.append(path[0])
        return [(path[i], path[i+1]) for i in range(self.word_length)]
    
    def list_symbols_strands(self):
        lstrands = self.list_strands()
        rd = self.rev_dict_lex_permutations()
        return [rd[item[0]][0]+rd[item[1]][0] for item in lstrands]

    def crossing(self, integer_pair):
        assert integer_pair in self.list_strands(), "pair not allowed"
        
    def get_list_numerical_pairs(self):
        dd = self.dict_lex_permutations()
        ll_num = [(dd[item[0]], dd[item[1]]) for item in self.list_LL]
        lr_num = [(dd[item[0]], dd[item[1]]) for item in self.list_LR]
        rl_num = [(dd[item[0]], dd[item[1]]) for item in self.list_RL]
        rr_num = [(dd[item[0]], dd[item[1]]) for item in self.list_RR]
        return ll_num, lr_num, rl_num, rr_num
    
    def sign_of_crossing(self):
        d = {}
        for first_pair in self.integer_pairs():
            for second_pair in self.integer_pairs():
                if (first_pair in self.ll_num and second_pair in self.rl_num):   
                    if end(second_pair)<end(first_pair):
                        d[(first_pair, second_pair)]=1
                    else: 
                        d[(first_pair, second_pair)]=0
                    continue
                    
                elif (first_pair in self.lr_num and second_pair in self.rl_num):   
                    d[(first_pair, second_pair)]=1
                    continue
                elif (first_pair in self.lr_num and second_pair in self.rr_num):     
                    if end(second_pair)<end(first_pair):
                        d[(first_pair, second_pair)]=1
                        continue
                    else: 
                        d[(first_pair, second_pair)]=0
                        continue
                elif (first_pair in self.rl_num and second_pair in self.lr_num):     
                    d[(first_pair, second_pair)]=-1
                    continue
                elif (first_pair in self.rl_num and second_pair in self.ll_num):     
                    if end(second_pair)>end(first_pair):
                        d[(first_pair, second_pair)]=-1
                        continue
                elif (first_pair in self.rr_num and second_pair in self.lr_num):     
                    if end(second_pair)>end(first_pair):
                        d[(first_pair, second_pair)]=-1
                        continue
                else:
                     d[(first_pair, second_pair)]=0
        return d
    
    def path_of_crossings(self):
        d = self.sign_of_crossing()
        l = self.list_strands()
        init = l[0]
        
    def symbol(self):
        d = {}
        for i, strand in enumerate(self.list_strands()):
            d[strand]=self.list_symbols_strands()[i]
        return d
                    
    def crossing_symbol(self, pair, symbol):
        crossings = []
        if symbol=='LL':
            for item in self.rl_num:
                if end(item)<end(pair):
                    crossings.append(1)
        if symbol=='LR':
            for item in self.rl_num:
                crossings.append(1)
            for item in self.rr_num:
                if end(item)<end(pair):
                    crossings.append(1)
        if symbol=='RL':
            for item in self.lr_num[::-1]:
                crossings.append(0)
            for item in self.ll_num[::-1]:
                if end(item)>end(pair):
                    crossings.append(0)
                
        if symbol=='RR':
            for item in self.lr_num[::-1]:
                if end(item)>end(pair):
                    crossings.append(0)
                
        return crossings
    
    
    def list_of_crossings_pairs(self):
        d = self.symbol()
        cross = self.sign_of_crossing()
        res = []
        for strand in self.list_strands():
            if d[strand]=='LL':
                for item in self.rl_num:
                    if (strand, item) in cross.keys() and cross[strand, item]!=0:
                        res.append((strand,item))
            if d[strand]=='LR':
                for item in self.rl_num:
                    if (strand, item) in cross.keys() and cross[strand, item]!=0:
                        res.append((strand,item))
                        
                for item in self.rr_num:
                    if (strand, item) in cross.keys() and cross[strand, item]!=0:
                        res.append((strand,item))
            if d[strand]=='RL':
                for item in self.lr_num[::-1]:
                    if (strand, item) in cross.keys() and cross[strand, item]!=0:
                        res.append((strand,item))
                for item in self.ll_num[::-1]:
                    if (strand, item) in cross.keys() and cross[strand, item]!=0:
                        res.append((strand,item))
            if d[strand]=='RR':
                for item in self.lr_num[::-1]:
                    if (strand, item) in cross.keys() and cross[strand, item]!=0:
                        res.append((strand,item))
                    
        return res

    def dowker_code(self):
        all_crossings = self.list_of_crossings_pairs()   
        list_sets = [sorted(item) for item in all_crossings]
        d =self.sign_of_crossing()
        cross = OrderedDict()
        for i,item in enumerate(all_crossings):
            cross[tuple(sorted(item))]=[]

        for i,item in enumerate(all_crossings):
            if (i+1)%2==0 and d[item]==1:
                cross[tuple(sorted(item))].append(-(i+1))
            elif (i+1)%2==0 and d[item]==-1:
                cross[tuple(sorted(item))].append(i+1)
            elif (i+1)%2 ==1:
                cross[tuple(sorted(item))].append(i+1)
                
        cv = cross.values()
        cv_sorted=[odd_order(item) for item in cv]    
        cvs = sorted(cv_sorted)
        dowker = [item[1] for item in cvs]
        
        return np.array(dowker)
    
    def convert_to_braid(self):
        self.ll_copy, self.lr_copy, self.rl_copy, self.rr_copy = self.get_list_numerical_pairs()
        self.remove_LR_RR_intersections()
        self.remove_LL_RL_intersections()
        #print self.ll_copy, self.lr_copy, self.rl_copy, self.rr_copy
        self.remove_LR_RL_intersections()
        #print self.ll_copy, self.lr_copy, self.rl_copy, self.rr_copy

        self.braid = self.permut_LR_RR+self.permut_LL_RL+self.permut_LR_RL
        return self
    
    def remove_LR_RR_intersections(self):
        #self.ll_copy, self.lr_copy, self.rl_copy, self.rr_copy = self.get_list_numerical_pairs()
        self.permut_LR_RR =[]
        if len(self.rr_copy)==0 or len(self.lr_copy)==0:
            return self
    
        self.max_lr = max([end(pair) for pair in self.lr_copy])
        self.min_rr = min([end(pair) for pair in self.rr_copy])
    
        while self.min_rr<self.max_lr:
            self.move_left_lr_rr()
    
        return self

    def move_left_lr_rr(self):
        for i,item_lr in enumerate(self.lr_copy):
            for j,item_rr in enumerate(self.rr_copy):
                if end(item_lr)==end(item_rr)+1:
                    self.lr_copy[i]=(start(item_lr), end(item_lr)-1)
                    self.rr_copy[j]=(start(item_rr), end(item_rr)+1)
                    self.permut_LR_RR.append(end(item_rr))
                    self.max_lr = max([end(pair) for pair in self.lr_copy])
                    self.min_rr = min([end(pair) for pair in self.rr_copy])
                    return self
        return self
    
    def remove_LL_RL_intersections(self):
        #self.ll_copy, self.lr_copy, self.rl_copy, self.rr_copy = self.get_list_numerical_pairs()
        self.permut_LL_RL =[]
        if len(self.ll_copy)==0 or len(self.rl_copy)==0:
            return self
    
        self.max_ll = max([end(pair) for pair in self.ll_copy])
        self.min_rl = min([end(pair) for pair in self.rl_copy])
    
        while self.min_rl<self.max_ll:
            self.move_right_ll_rl()
    
        return self

    def move_right_ll_rl(self):
        for i,item_ll in enumerate(self.ll_copy):
            for j,item_rl in enumerate(self.rl_copy):
                if end(item_ll)==end(item_rl)+1:
                    self.ll_copy[i]=(start(item_ll), end(item_ll)-1)
                    self.rl_copy[j]=(start(item_rl), end(item_rl)+1)
                    self.permut_LL_RL.append(end(item_rl))
                    self.max_ll = max([end(pair) for pair in self.ll_copy])
                    self.min_rl = min([end(pair) for pair in self.rl_copy])
                    return self
        return self
    
    def remove_LR_RL_intersections(self):
        #self.ll_copy, self.lr_copy, self.rl_copy, self.rr_copy = self.get_list_numerical_pairs()
        self.permut_LR_RL =[]
        if len(self.lr_copy)==0 or len(self.rl_copy)==0:
            return self
    
        self.max_lr = max([end(pair) for pair in self.lr_copy])
        self.min_rl = min([end(pair) for pair in self.rl_copy])
    
        while self.min_rl<self.max_lr:
            self.move_left_lr_rl()
    
        return self

    def move_left_lr_rl(self):
        for i,item_lr in enumerate(self.lr_copy):
            for j,item_rl in enumerate(self.rl_copy):
                if end(item_lr)==end(item_rl)+1:
                    self.lr_copy[i]=(start(item_lr), end(item_lr)-1)
                    self.rl_copy[j]=(start(item_rl), end(item_rl)+1)
                    self.permut_LR_RL.append(end(item_rl))
                    #if start(self.lr_copy[i])==end(self.lr_copy[i]):
                    #    self.lr_copy.pop(i)
                    #if len(self.lr_copy)==0:
                    #    self.max_lr=0
                    #    self.min_rl=1
                    #    return self
                    self.max_lr = max([end(pair) for pair in self.lr_copy])
                    self.min_rl = min([end(pair) for pair in self.rl_copy])
                    return self
        return self
    
    
if __name__ == '__main__':
    knot_string = 'LRLRRRLRRR'
    lk = LorenzKnot(knot_string)
    knot_code = str(lk.dowker_code())
    print "The Lorenz Knot with word ", knot_string, " has Dowker code \n",  knot_code 
    lk.convert_to_braid()
    print lk.braid