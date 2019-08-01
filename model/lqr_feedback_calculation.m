A = [[0 1];[0 0]]
B = [0;1]
Q = eye(2)
Q(1,1) = 1
R = 1

[K, S, CLP] = lqr(A,B,Q,R)