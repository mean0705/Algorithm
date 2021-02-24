from collections import Counter

def Change(denomination, total_coins):
    result = [total_coins+1] * (total_coins+1)
    coins_results = [[] for _ in range(total_coins+1)]
    result[0] = 0
    for i in range(1, total_coins+1):
        for coin in denomination:
            if i >= coin and result[i - coin] + 1 < result[i]:
                result[i] = result[i-coin] + 1
                coins_results[i] = coins_results[i-coin] + [coin]
    if result[total_coins] == total_coins+1:
        return []
    return coins_results[total_coins]

def PrintResult(result, index, denomination):
    print("Case", index)
    for d in denomination:
        if Counter(result)[d] != 0:
            print(d, "元", Counter(result)[d], "枚")
    print("\r")

def main():
    index = 1
    while 1:
        total_coins = int(input())
        if total_coins == 0:
            break
        denomination = [int(n) for n in input().split()]
        if index == 1:
            print("\n")
        PrintResult(sorted(Change(denomination, total_coins)), index, denomination)
        index += 1

if __name__ == "__main__":
    main()