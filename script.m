%% test1: calculate inverse of A
A = [1 2; 3 4]
Y = [1 0; 0 1]
X = GECP(A, Y)
abs(A*X - Y)
%% test2: calculate inverse of big A
n = 100;
A = rand([n,n]);
Y = eye(n);
X = GECP(A, Y);
max(abs(A*X - Y), [], "all")
%% test3: compare with linsolve
n = 100;
A = rand([n,n]);
Y = eye(n);
X1 = GECP(A, Y);
X2 = linsolve(A, Y);
max(abs(X1 - X2), [], "all")
%% test4: complex
n = 2;
A = rand([n,n]);
B = rand([n,n]);
a = rand([n,1]);
b = rand([n,1]);
C = A + 1i*B;
c = a + 1i*b;
M = [A, -B; B, A];
m = [a; b];
z2 = GECP(M, m);
z2 = z2(1:n) + 1i*z2(n+1:2*n);
z1 = linsolve(C, c);
max(abs(z1 - z2))
%% test5: complex GEPP
n = 2;
A = rand([n,n]);
B = rand([n,n]);
a = rand([n,1]);
b = rand([n,1]);
C = A + 1i*B;
c = a + 1i*b;
M = [A, -B; B, A];
m = [a; b];
z2 = GEPP(M, m);
z2 = z2(1:n) + 1i*z2(n+1:2*n);
z1 = linsolve(C, c);
max(abs(z1 - z2))
