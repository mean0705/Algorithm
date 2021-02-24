class Node:
    def __init__(self, freq, char):
        self.freq = freq
        self.char = char
        self.left = None
        self.right = None
        self.parent = None       

def createNodes(char, freq):
    newNode = Node(freq, char)
    return newNode

def createHuffmanTree(nodes): # huffman algorithm
    nodeQueue = nodes[:]
    while len(nodeQueue) > 1:
        nodeQueue.sort(key=lambda item:item.freq)
        node_left = nodeQueue.pop(0)
        node_right = nodeQueue.pop(0)
        node_parent = Node(node_left.freq + node_right.freq, None)
        node_parent.left = node_left
        node_parent.right = node_right
        node_left.parent = node_parent
        node_right.parent = node_parent
        nodeQueue.append(node_parent)
    nodeQueue[0].parent = None
    return nodeQueue[0]

def huffmanEncoding(nodes,root):
    codes = [''] * len(nodes)
    for i in range(len(nodes)):
        node_temp = nodes[i]
        while node_temp != root: 
            if node_temp.parent.left == node_temp:  codes[i] = '0' + codes[i]             
            else:                                   codes[i] = '1' + codes[i]                
            node_temp = node_temp.parent

    return codes

def huffmanDecoding(root, decodeString):
    answer = ""
    curr = root; 
    for i in range(len(decodeString)):
        if decodeString[i] == '0':  curr = curr.left           
        else:                       curr = curr.right
            
        if curr.left == None and curr.right == None: # is leaf, so we know which char
            answer += curr.char
            curr = root
    return answer

if __name__ == '__main__':
    case = 1
    print('Enter your input: ')
    while 1:
        num0fChar = int(input())
        if num0fChar == 0: break 
        listOfChar_Freq = []
        for _ in range(num0fChar):
            tempList = input().split()
            tempList[1] = int(tempList[1])
            listOfChar_Freq.append(tuple(tempList))

        decodeString = input()	
        nodes = []
        for j in range(len(listOfChar_Freq)):
            nodes.append(createNodes(listOfChar_Freq[j][0],listOfChar_Freq[j][1]))

        root = createHuffmanTree(nodes)
        codes = huffmanEncoding(nodes,root)
        final = []
        for k in zip(listOfChar_Freq,codes):
            final.append(k)
        final.sort(key=lambda x: x[0])


        print("\nHuffman Codes : #%d" %case)
        for item in final:
            print(item[0][0], item[1])
        print("Decode =", huffmanDecoding(root, decodeString))
        print("\r")
        case += 1