

def isEuler(visited, queue ):
    allVisited = True
    for e in visited:
        if e == 0:
            allVisited = False
    if allVisited:
        if queue[0] == queue[len(queue) - 1]:
            return 1
        else:
            return 2
    return 0

def printPath(flag, queue, result):
    resultStr = '\0'
    if flag == 1:
        for i in range(len(queue)):
            if i < len(queue) - 1:
                resultStr = resultStr + str(queue[i]) + " -> " 
            else:
                resultStr = resultStr + str(queue[i])
        result.append( resultStr )





# 搜尋過程只儲存一條路的狀態的資訊，搜尋結束後queue,visited會恢復為初始狀態
def dfs( u, visited, queue,  eulerEdges, result):
    queue.append(u)
    flag = isEuler( visited, queue )  # 判斷當前路徑是不是尤拉路，如果是則列印
    if flag > 0:
        printPath( flag, queue, result )
        
    for i in range(len(eulerEdges)):
        if visited[i] == 1:
            continue
        edge = eulerEdges[i]
        if edge[0] == u:
            visited[i] = 1
            dfs(edge[1], visited, queue,  eulerEdges, result )
            queue.pop()  # 將搜尋過的點彈出佇列
            visited[i] = 0  # 重置訪問狀態
        elif edge[1] == u:
            visited[i] = 1
            dfs(edge[0], visited, queue,  eulerEdges, result )
            queue.pop()  # 將搜尋過的點彈出佇列
            visited[i] = 0  # 重置訪問狀態

   
  
def main():
    case = 1
    print( "input:")
    eulerEdges = []
    start = 1  # 如果是尤拉道路必須從奇點開始
    queue = [] # 儲存路徑資訊
    result = []
    while 1:
        num, column = [int(n) for n in input().split()]
        if  num == 0 and column == 0:
            break
        else:
            eulerEdges.clear()
            result.clear()
            queue.clear()
            
            for i in range( column ):
                v1, v2 = [int(n) for n in input().split()]
                eulerEdges.append( [v1, v2] )
            
            visited = [0 for i in range(len(eulerEdges))] #訪問過的路
            dfs(start, visited, queue, eulerEdges, result )
            
            print( '\nCase:', case )
            if( len(result) > 0 ) :
                print( result[0] )
            else:
                print( 'No Euler Tours\n' )
            case += 1

if __name__ == '__main__':
    main()

 
'''
3 3
1 2 
1 3 
2 3 
6 10 
1 2 
1 3 
2 3 
2 4
2 5 
3 4
3 5 
4 5 
4 6 
5 6 
4 5 
1 2 
1 3 
1 4 
2 4 
3 4 
0 0
'''