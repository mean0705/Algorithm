
## 10627205_演算法第一次作業(4.數讀)

import numpy as np

def is_number(str):
  try:
    # 因為使用float有一個例外是'NaN'
    if str=='NaN':
      return False
    float(str)
    return True
  except ValueError:
    return False

class BSTNode:

	def __init__( self ):
		self.key   = None
		self.left  = None
		self.right = None


class BinarySearchTree:

	def __init__( self ):
		self.root = None
		self.n = 0

	def Insert( self, key ):
		ptr = BSTNode()
		ptr.key = key
		ptr.left = None
		ptr.right = None

		if ( self.root == None ):
			self.root = ptr
		else:
			parent = None
			current = self.root
			while ( current != None ):
				parent = current
				if ( key < current.key ):
					current = current.left
				else:
					current = current.right
			if ( key < parent.key ):
				parent.left = ptr
			else:
				parent.right = ptr
		self.n += 1

	def Postorder( self ):
		if ( self.root == None ):
			print( "No keys in the binary search tree." )
		else:
			print( "Binary Search Tree (Postorder):" )
			self._Postorder( self.root )
			print()

	def _Postorder( self, ptr ):
		if ( ptr != None ):
			self._Postorder( ptr.left )
			self._Postorder( ptr.right )
			print( ptr.key, end = " " )

	def Preorder( self ):
		if ( self.root == None ):
			print( "No keys in the binary search tree." )
		else:
			print( "Binary Search Tree (Preoder):" )
			self._Preorder( self.root )
			print()

	def _Preorder( self, ptr ):
		if ( ptr != None ):
			print( ptr.key, end = " " )
			self._Preorder( ptr.left )	
			self._Preorder( ptr.right )



if __name__ == "__main__" :

	while True :
	  node_list = list(input().split());
	  if node_list[0] == '0' and len(node_list) == 1 : 
		  break
	  tree = BinarySearchTree();
	  for j in range(0, len(node_list)) :
		  if (is_number(node_list[j])) : 
		      tree.Insert( int (node_list[j]) );
		  else :
			  tree.Insert( (node_list[j]) );
	  tree.Postorder();

	print('end');
