
import numpy as np

def f(num):
    sum = 0;
    for i in range(1,num+1) :
        sum = sum + i;
    return sum;


if __name__ == "__main__" :
    while True:
        num = int(input());
        if num == 0 : break;

        if num < 3 : print('error! we cannot find this case.');
        else :
            sum = 0;
            num = num - 3;
            while num > 0 :
                sum = sum + f(num);
                num = num - 2;
            print('this is answer', sum);

    print('end');
	

6
a 45
b 13
c 12
d 16
e 9
f 5
01001101
0

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
0 0
