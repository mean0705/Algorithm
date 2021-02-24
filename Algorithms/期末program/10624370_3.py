class Graph: 
  
    def __init__(self,vertices): 
        self.V= vertices #No. of vertices 
        self.graph = [] # default dictionary  
                                # to store graph 
          
   
    # function to add an edge to graph 
    def addEdge(self,u,v,w): 
        self.graph.append([u,v,w]) 
  
    # A utility function to find set of an element i 
    # (uses path compression technique) 
    def find(self, parent, i): 
        if parent[i] == i: 
            return i 
        return self.find(parent, parent[i]) 
  
    # A function that does union of two sets of x and y 
    # (uses union by rank) 
    def union(self, parent, rank, x, y): 
        xroot = self.find(parent, x) 
        yroot = self.find(parent, y) 
  
        # Attach smaller rank tree under root of  
        # high rank tree (Union by Rank) 
        if rank[xroot] < rank[yroot]: 
            parent[xroot] = yroot 
        elif rank[xroot] > rank[yroot]: 
            parent[yroot] = xroot 
  
        # If ranks are same, then make one as root  
        # and increment its rank by one 
        else : 
            parent[yroot] = xroot 
            rank[xroot] += 1
  
    # The main function to construct MST using Kruskal's  
        # algorithm 
    def KruskalMST(self, case): 
  
        result =[] #This will store the resultant MST 
        sum = 0
        i = 0 # An index variable, used for sorted edges 
        e = 0 # An index variable, used for result[] 
  
            # Step 1:  Sort all the edges in non-decreasing  
                # order of their 
                # weight.  If we are not allowed to change the  
                # given graph, we can create a copy of graph 
        self.graph =  sorted(self.graph,key=lambda item: item[2]) 
        parent = [] ; rank = [] 
  
        # Create V subsets with single elements 
        for node in range(self.V + 1 ):  #node start from 1
            parent.append(node) 
            rank.append(0) 
      
        # Number of edges to be taken is equal to V-1 
        while e < self.V -1 : 
  
            # Step 2: Pick the smallest edge and increment  
                    # the index for next iteration 
            u,v,w =  self.graph[i] 
            i = i + 1
            x = self.find(parent, u) 
            y = self.find(parent, v) 
  
            # If including this edge does't cause cycle,  
                        # include it in result and increment the index 
                        # of result for next edge 
            if x != y: 
                e = e + 1     
                result.append([u,v,w]) 
                self.union(parent, rank, x, y)             
            # Else discard the edge 
  
        # print the contents of result[] to display the built MST 
        for u,v,weight  in result: 
            sum = sum + weight
        print( '\n\nCase' , case ,':') 
        print( 'Minimum Cost:', sum ,'\n' )


if __name__ == '__main__':
    case = 1
    print( "input:")
    while 1:
	    line = input().split()
		num = int(line[0])
		column = int(line[1])
        if  num == 0 and column == 0:
            break
        else:
            g = Graph( num )
            for i in range( column ):
                line = input().split()
		        v1 = int(line[0])
		        v2 = int(line[1])
				weight = int(line[2])
                g.addEdge(v1, v2, weight)
            g.KruskalMST( case )
            case += 1
    
'''
4 4
1 2 10 
1 3 8
2 4 5 
3 4 2 
5 7 
1 2 2 
1 4 10 
1 5 6 
2 3 5 
2 5 9 
3 5 8 
4 5 12 
00
'''
4 4
1 2 10
1 3 8
2 4 5
3 4 2
5 7
1 2 2
1 4 10
1 5 6 
2 3 5
2 5 9
3 5 8
4 5 12
00