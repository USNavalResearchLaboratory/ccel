#cython: embedsignature=True
# Author: Michael Shlanta
# michael.shlanta@nrl.navy.mil, michael.shlanta@dsu.edu
# This is the Cython implementation of the Corrected Conditional Entropy (CCE) detector
# which binds the C implementation into a Python library 
from libCCE cimport llRoot, createLinkedList, insertSequence, freeLL, tNode, createTree, printLayer, printList, calcCCEs, applyToList, llNode, isnan, resetArena
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.stdlib cimport free, malloc
from libc.stdio cimport printf
from collections import Counter
from scipy.stats import entropy
from itertools import islice
import numpy as np
import cython

@cython.embedsignature(True)
cdef class CCE(object):
    """
    This class is a wrapper for the C CCE implementation which is based
    on the work done by S. Gianvecchio and H. Wang's article on detecting
    covert timing channels using an entropy-based approach. In particular,
    please refer to section 3.4 (this section describes CCE implementation using
    a q-ary tree of height m) in the following publication: 
    S. Gianvecchio, H. Wang, "Detecting Covert Timing Channels: An entropy-based
    Approach", Proceedings of CCS'07, Oct.29-Nov.2, 2007.  
    
    This class calculates Corrected Conditional Entropy (CCE) of a given sequence.  
    The inital body of samples can be provided to `__init__` and they 
    will be populated into the tree using a faster pure C method.  
    Any new data can be provided to `insertSequence` which will insert 
    it into the tree.

    As a general rule this implementation will try and fail via 
    assertion before returining compromised information to the user.
    """
    # instance level c attrubutes 
    cdef tRoot* root
    cdef unsigned int branchingFactor
    
    def __cinit__(self, unsigned int branchingFactor=5, initalSequences=None, int subSeqLen=50):
	"""
	Called on creation of object, if initalSequences is supplied, the sequences will be 
  	populated to the tree via a faster version of insertSequence.

	Args:
        	branchingFactor (int): 
            		The branching factor of the tree, this is equivalent to
            		the number of bins you have.
        	initalSequences (iterable):
            		A python iterable which contains Inter-Packet Delays(IPDs) for insertion.
        	subSeqLen (int): 
            		The length you would like the subsequences devided into.
    	Raises:
        	ValueError:
            		if subSeqLen is less than the length of the original 
            		sequence or if the branching factor is less than or equal to zero.
	"""
        cdef int seqLen 
        cdef int* cSeq
        cdef int i = 0
        cdef int offset

        if branchingFactor <= 0:
            raise ValueError('branching factor must be > 0.')

        self.root = createTree(branchingFactor)
        self.branchingFactor = branchingFactor

        if initalSequences is not None:
            # allocate memory to pass the values from the python object to the CCE code
            seqLen = len(initalSequences)
            if subSeqLen > seqLen:
                raise ValueError("subSeqLen must be less than the length of the initial sequance")
            cSeq = <int*>PyMem_Malloc(seqLen * sizeof(int))
            
            # Copy the python object to the memory
            for seq in initalSequences:
                cSeq[i] = seq
                if cSeq[i] >= branchingFactor:
                    raise ValueError("Bin value %d at index %d is not in range [0, %d)" % (seq, i, self.branchingFactor))
                i += 1
            
            # insert sequences into the tree using pointer math to
            # set the section we are interested in.
            for i in range((seqLen - subSeqLen)):
                insertSequence(self.root, cSeq + i, subSeqLen)
            
            PyMem_Free(cSeq)


    def __init__(self, branchingFactor=5, initalSequence=None, subSeqLen=50):
        """
	All real work is done in __cinit__
	"""
        pass

    def __dealloc__(self):
        """ 
        Used to cleanup memory on destruction of object. Automatically
        called. Should NOT be manually called it as it WILL break everything.
        """
        freeTree(self.root)

    def insertSequence(self, pySeq):
        """
        Wrapper for insertSequence, converts a supplied python
	iterable which supports indexing and len() into an integer
	buffer of the same length.

        Args:
            pySeq (list or tuple):
                A python iterable that supports indexing and len()
        """
        cdef int* cSeq
        cdef int i
        cdef int seqLen

        seqLen = len(pySeq)
        cSeq = <int*>PyMem_Malloc(seqLen * sizeof(int))
        # populate the C buffer with values from the python object.
        for i in range(seqLen):
            cSeq[i] = pySeq[i]
            if cSeq[i] >= self.branchingFactor:
                raise ValueError("Bin value %d at index %d is not in range [0, %d)" % (cSeq[i], i, self.branchingFactor))
        insertSequence(self.root, cSeq, seqLen)
        PyMem_Free(cSeq)

    def calculateCCE(self):
        """
        Calculate the CCE for all inserted sequences and return the
        CCE for all sub-sequence lengths upto the maximum sub-
        sequence length.

        Returns:
            CCEs(list): 
		A list of floats containing the CCE for each layer
		of the tree. This list will only be as long as the
		number of unique sequences.
        """
        cdef double* CCEArr
        cdef llNode* n
        cdef int i = 0

        CCEs = []

        CCEArr = calcCCEs(self.root)
        
        # The array returned from calcCCEs is NaN terminated so we
        # loop till we find the NaN.
        while not isnan(CCEArr[i]):
            CCEs.append(CCEArr[i])
            i += 1

        # free memory allocated in calcCCEs to hold the return values.
        free(CCEArr)

        return CCEs

    def resetTree(self):
        '''
        This function can be used to reset the tree and reuse it 
        by freeing all the memory assocated with the tree.

        This can be very useful when you want to calculate new data
        that is of a similar size to the previously-analized data.
        '''
        resetArena(self.root.arenaManagement)
#end of class CCE definition

def _window(seq, int n=2):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

cdef double _precentUnique(frequencys):
    f = [x for x in frequencys if x == 1] 
    return <double>(len(f) / sum(frequencys))

@cython.embedsignature(True)
def localMinCCE(seq, int maxLen):
    '''
    This is a faster (although most likely less accurate) implementation of CCE,
    which only computes a local minimum. 

    It is based on the following publication:
    A. Porta, G. Baselli, D. Liberati, N. Montano, C. Cogliati, T.
    Gnecchi-Ruscone, A. Malliani, and S. Cerutti, “Measuring
    Regularity by Means of a Corrected Conditional Entropy in
    Sympathetic Outflow,” Biological Cybernetics, vol. 78, no. 1, 
    pp. 71-78, Jan. 1998. 
    
    This function calculates the CCE of a sequence using the local minimum as an
    early stopping metric. It still uses the same C library and the
    the same tree approach to calculating CCE, but stops when it finds
    a local minimum. 

    Args:
        seq (iterable):
            An iterable containing the pre-binned values of your IPD
            sequence.
        maxLen (int):
            The maximum subsequence length to be tried.  This is
            more of an upper bound as the local minimum is usually
            found before this is reached.
    
    Returns:
        minCCE (float):
		The local minimum CCE of the provided sequence.
    '''
    cdef int seqLen 
    cdef double firstOrderEntropy
    cdef double correctionFactor
    cdef double previousEntropy
    cdef double uniques
    cdef double minCCE
    cdef double CCE
    cdef double ent

    previousEntropy = 0.0
    firstOrderEntropy = 0.0
    minCCE = float('inf')
    # print('CCE, Conditional Entropy, Correction Factor') #debug info
    
    # Prime the first order entropy.
    c = Counter(_window(seq[:len(seq)-maxLen], 1))
    frequency = list(c.values())
    ent = entropy(frequency)
    firstOrderEntropy = ent

    for seqLen in range(2, maxLen+1):
        conditionalEntroy = ent - previousEntropy
        previousEntropy = ent
        uniques = _precentUnique(frequency)
        correctionFactor = firstOrderEntropy * uniques
        CCE = conditionalEntroy + correctionFactor
        # print("%f, %f, %f" % (CCE, conditionalEntroy, correctionFactor)) # debug info
        if CCE > minCCE:
            break
        else:
            minCCE = CCE
        if uniques == 1.0:
            break
        # This kind of slicing ensures that the output of the tree-based
        # version and this one match.  You lose out a few
        # samples per length, but they would have been from an
        # incomplete segment anyway.
        c = Counter(_window(seq[:len(seq)-(maxLen-(seqLen-1))], seqLen))
        frequency = list(c.values())
        # print(frequency) # debug info
        ent = entropy(frequency)

    return minCCE
