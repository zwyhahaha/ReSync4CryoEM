function Y = solve_SDP(S)
[n, ~] = size(S);
n = n / 2;
cvx_begin sdp
variable X(2*n,2*n) symmetric;
maximize trace(S*X) 
subject to 
X>=0;%==semidefinite(2*n);
diag(X)==ones(2*n,1);
diag(X(n+1:2*n, 1:n))==zeros(n,1);
diag(X(1:n, n+1:2*n))==zeros(n,1);
cvx_end 
Y = X;
end 