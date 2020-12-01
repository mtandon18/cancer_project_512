import math

def calc_modifier(logFChange, i, j):
    M = abs((logFChange[i] + logFChange[j])/2.0)
    N = abs(logFChange[i] - logFChange[j])
    return 1 - 2.0/math.pi*math.atan(M/3.0) + math.asinh(N/3.0)/2.0