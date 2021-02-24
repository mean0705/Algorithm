#
#	AlgorithmsTest.py 
#	Copyright@ Yuan-Hsiang Chang, Ph.D.
#	Department of Information & Computer Engineering
#	Chung Yuan Christian University
#
from Algorithms import *

import numpy as np
import random
import time

#
#  排列組合
#
A = [ 1, 2, 3, 4 ]    # 或 A = np.array( [ 1, 2, 3, 4 ] )
print( "Subsets of A" )
Subsets( A )

print( "Permutations of A" )
Permutations( A )


#
#  排序演算法
#
A1 = [ i + 1 for i in range( 10000 ) ]
random.shuffle( A1 )

A = A1.copy()
start_time = time.time()
BubbleSort( A )
print( "Bubble Sort --- %s seconds ---" % ( time.time() - start_time ) )

A = A1.copy()
start_time = time.time()
MergeSort( A )
print( "Merge Sort --- %s seconds ---" % ( time.time() - start_time ) )

A = A1.copy()
start_time = time.time()
QuickSort( A )
print( "Quick Sort --- %s seconds ---" % ( time.time() - start_time ) )


#
#  圖形演算法
#
G = AdjacencyMatrix( 6 )
G.SetEdge( 1, 2 )
G.SetEdge( 1, 4 )
G.SetEdge( 1, 5 )
G.SetEdge( 2, 3 )
G.SetEdge( 2, 5 )
G.SetEdge( 3, 6 )
G.SetEdge( 5, 6 )
G.Display()
G.BFS( 1 )
G.DFS( 1 )