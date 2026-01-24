function [time,err] = TimeMethods(A,B,a,b)
%TIMEMETHODS Tests different Gaussian Elimination methods on complex
%equation system
%   Test on complex equation system Cz = c
%   where C = A + Bi, c = a + bi
arguments (Input)
    A % nxn matrix (real part)
    B % nxn matrix (imaginary part)
    a % nx1 vector (real part)
    b % nx1 vector (imaginary part)
end

arguments (Output)
    time % 5x1 vector, elapsed time for each method
    err % 5x1 vector, error for each method
end

n = size(A, 1);
C = A + 1i*B;
c = a + 1i*b;
M = [A,-B;B,A];
m = [a; b];
time = zeros(5,1);
err = zeros(5,1);
p = 2;

tic
z = linsolve(C, c); % builtin solver
time(1) = toc;
err(1) = norm(C*z-c, p);

tic
z = GEPP(C,c); % complex GEPP
time(2) = toc;
err(2) = norm(C*z-c, p);

tic
z = GECP(C,c); % complex GECP
time(3) = toc;
err(3) = norm(C*z-c, p);

tic
z = GEPP(M, m); % real block GEPP
z = z(1:n) + 1i*z(n+1:2*n);
time(4) = toc;
err(4) = norm(C*z-c, p);

tic
z = GECP(M, m); % real block GECP
z = z(1:n) + 1i*z(n+1:2*n);
time(5) = toc;
err(5) = norm(C*z-c, p);

end