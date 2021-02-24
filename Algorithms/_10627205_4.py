## 10627205_演算法第一次作業(4.數讀)

import numpy as np

def firstpos(board) :
    for x in range(9):
        for y in range(9) :
            if board[x][y] == 0 :
                return x, y
    return False, False

def nextpos(board, x, y) :
    for y2 in range(y+1, 9) :
        if board[x][y2] == 0:
            return x, y2
    for x2 in range(x+1, 9) :
        for y2 in range(0, 9) :
            if board[x2][y2] == 0 :
                return x2, y2
    return -1, -1

def findvalue(board, x, y) :
    i, j = x//3, y//3
    ninegrid = [board[i*3+n][j*3+m] for n in range(3) for m in range(3)]                           #當前九宮格
    vals = set([x for x in range(1, 10)]) - set(ninegrid) - set(board[x]) - set(list(zip(*board))[y])
    return list(vals)

def RunSudoku(board, x, y) :

    for vals in findvalue(board, x, y) :
        board[x][y] = vals
        x2, y2 = nextpos(board, x, y)
        if y2 == -1 :
            return True
        else :
            v2 = RunSudoku(board, x2, y2)
            if v2 :
                return True
            board[x][y] = 0                  #若沒有解開則重新找



if __name__ == "__main__" :

    board = [[0 for i in range(9)] for j in range(9)]
    for i in range(0,9) :
        print("請輸入第", i+1, "列的數字(9筆): " )
        da = list(input())
        board[i] = da
        for j in range(0, 9) :
            board[i][j] = int(board[i][j])
        del da
    x, y = firstpos(board)
    RunSudoku(board, x, y)
    print(board)