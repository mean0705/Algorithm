gold = 0

def findP(map) :
    for i in range(len(map)) :
        for j in range(len(map[i])) :
            if map[i][j] == 'P' :
                return i, j

def checkDirect(map, x, y) :
    if map[x][y+1] == '-' :
        #"Right"
        return x, y+1
    elif map[x+1][y] == '-' :
        #"Down"
        return x+1, y
    elif map[x][y-1] == '-' :
        #"Left"
        return x, y-1
    elif map[x-1][y] == '-' :
        #"Up"
        return x-1, y
    else :
        return -5, -5

def checkEvent(map, usrMap, x, y) :
    if map[x][y] == 'G' :
        global gold 
        gold += 1
    if map[x][y] == '#' :
        return "wall"
    elif map[x-1][y] == 'T' or map[x+1][y] =='T' or map[x][y+1] =='T' or map[x][y-1] =='T' :
        return "trap"
    else :
        return "normal"
    

def goldGame(map, pX, pY) :
    # create a empty map for user
    usrMap = [[] for i in range(len(map))]
    for i in range(len(map)) :
        for j in range(len(map[i])) :
            usrMap[i].append('-')
    # created
    usrMap[pX][pY] = 'P'

    currentX, currentY = pX, pY
    finish = False 
    while( finish == False ) :
        newX, newY = checkDirect(usrMap, currentX, currentY) 
        # check is it end
        if newX == -5 :
            if currentX != pX or currentY != pY :
                currentX, currentY = pX, pY
                continue
            elif currentX == pX and currentY == pY :
                finish = True
                continue
        # done
        event = checkEvent(map, usrMap, newX, newY)
        if event == 'trap' :
            temp = usrMap[currentX][currentY]
            usrMap[newX+1][newY] = 'T'
            usrMap[newX-1][newY] = 'T'
            usrMap[newX][newY-1] = 'T'
            usrMap[newX][newY+1] = 'T'
            usrMap[newX][newY] ='.'
            usrMap[currentX][currentY] = temp
        elif event == 'wall' :
            usrMap[newX][newY] = '#'
        elif event == 'normal' :
            usrMap[newX][newY] = '.'
            currentX, currentY = newX, newY
        else :
            print("Some thing wrong",event)
        
def main() :
    print("Start your input:")
    global gold
    while 1  :
        gold= 0
        length, width = [int(n) for n in input().split()]
        if length == 0 and width == 0 :
            break

        read = []

        for i in range( width ) :
            read.append(input())

        map = [[] for i in range(width) ] 
        for i in range( width ) :
            for item in str(read[i]) :
                map[i].append(item)
         # read input end
        pIndexX, pIndexY = findP(map)
        goldGame(map, pIndexX, pIndexY) 
        print("\n",gold)

if __name__ == '__main__':
    main()

'''
7 4 
#######
#P.GTG#
#..TGG#
#######
8 6 
######## 
#...GTG# 
#..PG.G# 
#...G#G# 
#..TG.G# 
######## 
0 0
'''
