class Node:
    def __init__(self, freq, char):
        self.left = None
        self.right = None
        self.father = None
        self.freq = freq
        self.char = char
    def isLeft(self):
        return self.father.left == self

def createNodes(freqs, char):
    newNode = Node(freqs, char)
    return newNode

def createHuffmanTree(nodes):
    queue = nodes[:]
    while len(queue) > 1:
        queue.sort(key=lambda item:item.freq)
        node_left = queue.pop(0)
        node_right = queue.pop(0)
        node_father = Node(node_left.freq + node_right.freq, None)
        node_father.left = node_left
        node_father.right = node_right
        node_left.father = node_father
        node_right.father = node_father
        queue.append(node_father)
    queue[0].father = None
    return queue[0]

def huffmanEncoding(nodes,root):
    codes = [''] * len(nodes)
    for i in range(len(nodes)):
        node_tmp = nodes[i]
        while node_tmp != root:
            if node_tmp.isLeft():
                codes[i] = '0' + codes[i]
            else:
                codes[i] = '1' + codes[i]
            node_tmp = node_tmp.father
    return codes

def huffmanDecoding(root, decodeString):
    answer = ""
    curr = root; 
    for i in range(len(decodeString)):
        if decodeString[i] == '0':
            curr = curr.left
        else:
            curr = curr.right
        if curr.left == None and curr.right == None:
            answer += curr.char
            curr = root
    return answer

if __name__ == '__main__':
    index = 1
    while 1:
        char_n = int(input())
        if char_n == 0:
            break
        chars_freqs = []
        for _ in range(char_n):
            tempList = input().split()
            tempList[1] = int(tempList[1])
            chars_freqs.append(tuple(tempList))
        decodeString = input()
        nodes = []
        for i in range(len(chars_freqs)):
            nodes.append(createNodes(chars_freqs[i][1],chars_freqs[i][0]))
        root = createHuffmanTree(nodes)
        codes = huffmanEncoding(nodes,root)
        temp = []
        for item in zip(chars_freqs,codes):
            temp.append(item)
        temp.sort(key=lambda tup: tup[0])
        if index == 1:
            print("\n")
        print("Huffman Codes #%d" %index)
        for item in temp:
            print(item[0][0], item[1])
        print("Decode =", huffmanDecoding(root, decodeString))
        print("\r")
        index += 1

'''
6
a 45
b 13
c 12
d 16
e 9
f 5
01001101
6
A 2
B 6
C 15
D 12
E 8
F 3
010101001100
0
'''