import numpy as np

def minDepth (list, zero_sp) :
    list1 = list[:];
    list2 = list[:];
    print(list);
    rDepth = 99;
    dDepth = 99;
    if zero_sp == 8 :
        if int(list[7]) == 1 :
            return int(list[0]) + int(list[1]) + int(list[3]) + int(list[4]) + int(list[5]) + int(list[6]) + 3;
        elif int(list[5]) == 1 :
            return int(list[0]) + int(list[1]) + int(list[2]) + int(list[3]) + int(list[4]) + int(list[7]) + 3;
    elif zero_sp == 2 or zero_sp == 5 :
        x = int(list[zero_sp + 3]);
        list1[zero_sp], list1[zero_sp + 3] = list1[zero_sp + 3], list1[zero_sp];
        dDepth = minDepth(list1, zero_sp +3) + x;
    elif zero_sp == 6 or zero_sp == 7 :
        y = int(list[zero_sp + 1]);
        list2[zero_sp], list2[zero_sp + 1] = list2[zero_sp + 1], list2[zero_sp];
        rDepth = minDepth(list2, zero_sp +1) + y;
    else :
        y = int(list[zero_sp + 1]);
        list2[zero_sp], list2[zero_sp + 1] = list2[zero_sp + 1], list2[zero_sp];
        rDepth = minDepth(list2, zero_sp +1) + y;

        x = int(list[zero_sp + 3]);
        list1[zero_sp], list1[zero_sp + 3] = list1[zero_sp + 3], list1[zero_sp];
        dDepth = minDepth(list1, zero_sp +3) + x;

    if rDepth < dDepth : return rDepth;
    else               : return dDepth;


if __name__ == "__main__" :
    mm = list(map(int, input().split()))
    nn = list(map(int, input().split()))
    xx = list(map(int, input().split()))
    my_diagram = mm + nn + xx;
    zero_sp = 0;
    ans = minDepth(my_diagram, zero_sp);
    print (ans);