function fi=UGEstatics(fi, dop, n , p, n_max, dx, delta_acc, x)


dx2 = dx*dx;

    a = ones(1,n_max)/dx2;
    c = ones(1,n_max)/dx2;
    b = -(2/dx2 + n + p);
    f = n - p - dop - fi.*(n + p);
    
%(B) Define the elements of the coefficient matrix and initialize the forcing
%    function at the ohmic contacts 

a(1) = 0;
c(1) = 0;
b(1) = 1;
f(1) = fi(1);
a(n_max) = 0;
c(n_max) = 0;
b(n_max) = 1;
f(n_max) = fi(n_max);

% define LU matrices for Ax=b where A=LU and find alpha, beta
%(C)  Start the iterative procedure for the solution of the linearized Poisson
%     equation using LU decomposition method:

flag_conv = 0;		           % convergence of the Poisson loop
k_iter= 0;
tic
while(~flag_conv)            
    k_iter = k_iter + 1 
figure (100); plot(x, fi); xlabel ('distance'); ylabel('normalized energy'); hold on;     
pause;

fiold=fi;  % storing previous fi to compare error
fi=LUdecomp3(a,b,c,f, n_max); % calling LU decomposition function

delta = fi - fiold;
    %delta_max=max(abs(delta))
    delta_max (k_iter) = max(abs(delta));
    
    
sprintf('delta_max = %d',delta_max(k_iter));      %'k_iter = %d',k_iter,'
    
    if(delta_max (k_iter) < 0.01*delta_acc)
        flag_conv = 1
    else
        % update functions b and f
        
%     b = -(2/dx2 + exp(fi) + exp(-fi));
%     f = exp(fi) - exp(-fi) - dop - fi.*(exp(fi) + exp(-fi));
        
p= exp(-fi); 
n= exp(fi);
        b = -(2/dx2 + n + p);
        f = n - p - dop - fi.*(n + p);
            % initialize boundaries
            a(1) = 0;
            c(1) = 0;
            b(1) = 1;
            f(1) = fi(1);
            a(n_max) = 0;
            c(n_max) = 0;
            b(n_max) = 1;
            f(n_max) = fi(n_max);
        
        
    end
  
%whos(k_iter, delta_max)
figure (101); semilogy(1:k_iter, delta_max, 'ro-');  xlabel ('interation'); ylabel('max error');hold on;

end
figure (100); hold off;  
