class Graph: 
  
    def __init__(self,vertices): 
        self.V = vertices
        self.graph = [] 

    def addEdge(self,person1,person2,weight): 
        self.graph.append([person1,person2,weight]) 

    def find(self, parent, i): 
        if parent[i] != i: return self.find(parent, parent[i]) 
        else: return i 
        
    def union(self, parent, rank, x, y): 
        xroot = self.find(parent, x) 
        yroot = self.find(parent, y) 
  
        if rank[xroot] < rank[yroot]: 
            parent[xroot] = yroot 
        elif rank[xroot] > rank[yroot]: 
            parent[yroot] = xroot 
        else : 
            parent[yroot] = xroot 
            rank[xroot] += 1

    def KruskalMST(self, case):  
        result =[]  
        parent = [] 
        rank = []
        self.graph = sorted(self.graph,key=lambda h: h[2]) 
  
        for node in range(self.V + 1 ): # start from 1
            parent.append(node) 
            rank.append(0) # number of connect
 
        v = 0
        i = 0 
        while v < self.V -1 : 
            person1,person2,weight =  self.graph[i] 
            i = i + 1
            x = self.find(parent, person1) 
            y = self.find(parent ,person2) 
  
            if x != y: 
                v = v + 1     
                result.append([person1,person2,weight]) 
                self.union(parent, rank, x, y)             

        sum = 0
        for person1,person2,weight in result: 
            sum = sum + weight
        print( '\n\nCase' , case ,':') 
        print( 'Minimum Cost:', sum ,'\n' )


if __name__ == '__main__':
    case = 1
    print('Enter your input: ')
    while 1:      
        line = input().split()
        numOfPerson = int(line[0])
        numOfData = int(line[1])
        if  numOfPerson == 0 and numOfData == 0:
            break
        else:
            g = Graph( numOfPerson )
            for i in range( numOfData ):
                line = input().split()
                person1 = int(line[0])
                person2 = int(line[1])
                weight = int(line[2])
                g.addEdge(person1, person2, weight)

            g.KruskalMST( case )
            case += 1
