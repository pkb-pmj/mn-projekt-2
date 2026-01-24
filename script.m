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

%%
n = 150;
iters = 150;
T = zeros([iters,5]);
E = zeros([iters,5]);

for i = 1:iters
    A = rand([n,n]);
    B = rand([n,n]);
    a = rand([n,1]);
    b = rand([n,1]);
    C = A + 1i*B;
    c = a + 1i*b;
    M = [A,-B;B,A];
    m = [a; b];
    t = zeros(5,1);
    e = zeros(5,1);
    p = 2;
    
    tic
    z = linsolve(C, c); % wbudowana funkcja
    t(1) = toc;
    e(1) = norm(C*z-c, p);
    
    tic
    z = GEPP(C,c); % GEPP na zespolonej
    t(2) = toc;
    e(2) = norm(C*z-c, p);
    
    tic
    z = GECP(C,c); % GECP na zespolonej
    t(3) = toc;
    e(3) = norm(C*z-c, p);
    
    tic
    z = GEPP(M, m); % GEPP z użyciem macierzy blokowej
    z = z(1:n) + 1i*z(n+1:2*n);
    t(4) = toc;
    e(4) = norm(C*z-c, p);
    
    tic
    z = GECP(M, m); % GECP z użyciem macierzy blokowej
    z = z(1:n) + 1i*z(n+1:2*n);
    t(5) = toc;
    e(5) = norm(C*z-c, p);

    T(i, :) = t;
    E(i, :) = e;
end

methodNames = ["linsolve","GEPP(C,c)","GECP(C,c)","GEPP(M, m)","GECP(M, m)"];

% plot time
% speed gradually increases during first 10-30 iterations, possibly JIT compilation?
%plot(1:iters, T)
%set(gca,'yscale','log')
%xlabel("iteration")
%ylabel("time")
%legend(methodNames)

% plot error
%plot(1:iters, E)
%set(gca,'yscale','log')
%xlabel("iteration")
%ylabel("error")
%legend(methodNames)

% plot error vs time
scatter(T, E)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("time")
ylabel("error")
xlim([1e-6,1e0])
ylim([1e-16 1e-12])
%legend(methodNames)


set(gcf, 'Position', [100 100 300 300])

% zapis do pliku PNG
filename = sprintf('error_time_n=%d.png', n);
print(gcf, filename, '-dpng', '-r300')   % 300 DPI
%% plot error correlation
% errors for all methods are closely correlated
plotmatrix(E);
ax = findall(gcf,'Type','axes');
% set(ax,'XScale','log','YScale','log')
for k = 1:numel(ax)
    if any(strcmp({ax(k).Children.Type},'line'))
        set(ax(k),'XScale','log','YScale','log')
    end
end

%% plot time correlation
% time for all methods is also correlated
% makes sense, if GEPP swaps rows then GECP also swaps rows
plotmatrix(T);
ax = findall(gcf,'Type','axes');
% set(ax,'XScale','log','YScale','log')
for k = 1:numel(ax)
    if any(strcmp({ax(k).Children.Type},'line'))
        set(ax(k),'XScale','log','YScale','log')
    end
end


%%
%% mean time vs matrix size n

ns = [20 40 60 80 100 120 140 160 200];
iters = 10;          % uśrednianie
meanT = zeros(length(ns), 5);
meanE = zeros(length(ns), 5);

for k = 1:length(ns)
    n = ns(k);
    T = zeros(iters,5);
    E = zeros(iters,5);
    for i = 1:iters
        A = rand([n,n]);
        B = rand([n,n]);
        a = rand([n,1]);
        b = rand([n,1]);
        C = A + 1i*B;
        c = a + 1i*b;
        M = [A,-B;B,A];
        m = [a; b];

        tic
        z = linsolve(C, c);
        T(i,1) = toc;
        E(i,1) = norm(C*z-c, p);

        tic
        z = GEPP(C,c);
        T(i,2) = toc;
        E(i,2) = norm(C*z-c, p);

        tic
        z = GECP(C,c);
        T(i,3) = toc;
        E(i,3) = norm(C*z-c, p);

        tic
        z = GEPP(M, m);
        z = z(1:n) + 1i*z(n+1:2*n);
        T(i,4) = toc;
        E(i,4) = norm(C*z-c, p);

        tic
        z = GECP(M, m);
        z = z(1:n) + 1i*z(n+1:2*n);
        T(i,5) = toc;
        E(i,5) = norm(C*z-c, p);
    end

    meanT(k,:) = mean(T,1);
    meanE(k,:) = mean(E,1);
end

%% plot mean time & error vs n
figure
plot(ns, meanT, '-o', 'LineWidth', 1.5)
set(gca,'yscale','log')
xlabel('rozmiar macierzy: n')
ylabel('średni czas[s]')
legend(methodNames, 'Location','northwest')
title('Zależność czasu działania metody od rozmiaru macierzy')
grid on

set(gcf, 'Position', [100 100 600 400])

% zapis do pliku PNG
filename = 'time_vs_n';
print(gcf, filename, '-dpng', '-r300')   % 300 DPI



figure
plot(ns, meanE, '-o', 'LineWidth', 1.5)
set(gca,'yscale','log')
xlabel('rozmiar macierzy: n')
ylabel('średni błąd')
legend(methodNames, 'Location','northwest')
title('Zależność błędu metody od rozmiaru macierzy')
grid on

set(gcf, 'Position', [100 100 600 400])

% zapis do pliku PNG
filename = 'error_vs_n';
print(gcf, filename, '-dpng', '-r300')   % 300 DPI
