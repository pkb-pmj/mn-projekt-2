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
%% test6: Czy trik ze macierzą blokową coś daje?
n = 3;
A = rand([n,n]);
B = rand([n,n]);
a = rand([n,1]);
b = rand([n,1]);
C = A + 1i*B;
c = a + 1i*b;
M = [A,-B;B,A];
m = [a; b];
t = zeros(5,1);
tic
z0 = linsolve(C, c); % wbudowana funkcja
t(1) = toc;
tic
z1 = GECP(C,c); % GECP na zespolonej
t(2) = toc;
tic
z2 = GECP(M, m); % GECP z użyciem macierzy blokowej
z2 = z2(1:n) + 1i*z2(n+1:2*n);
t(3) = toc;
tic
z3 = GEPP(C,c); % GEPP na zespolonej
t(4) = toc;
tic
z4 = GEPP(M, m); % GEPP z użyciem macierzy blokowej
z4 = z4(1:n) + 1i*z4(n+1:2*n);
t(5) = toc;

E = zeros(5,1);
p = Inf;
E(1) = norm(C*z0-c,p);
E(2) = norm(C*z1-c,p);
E(3) = norm(C*z2-c,p);
E(4) = norm(C*z3-c,p);
E(5) = norm(C*z4-c,p);

format short e
disp('Porównanie błędów ||C*z-c|| dla różnych metod')
disp(table( ...
    categorical(["linsolve";"GECP(C,c)";"GECP(M, m)";"GEPP(C,c)";"GEPP(M, m)"]), ...
    t, ...
    E, ...
    'VariableNames',{'Metoda','Czas',sprintf('norm(C*z-c,%g)',p)}))