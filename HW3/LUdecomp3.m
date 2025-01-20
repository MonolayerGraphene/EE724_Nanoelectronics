function fi=LUdecomp3(a,b,c,f, n_max)
alpha(1)=b(1);
for i=2:n_max
beta(i)=a(i)./alpha(i-1);
alpha(i)=b(i) - beta(i).*c(i-1);
end

% Solution of Lv = f %    

    v(1) = f(1);
    for i = 2:n_max
        v(i) = f(i) - beta(i)*v(i-1);
    end

    % Solution of U*fi = v %    

    temp = v(n_max)/alpha(n_max);
    %delta(n_max) = temp - fi(n_max);
     fi(n_max)=temp;
    for i = (n_max-1):-1:1       %delta%
        fi(i) = (v(i)-c(i)*fi(i+1))/alpha(i);
        %delta(i) = temp - fi(i);
        %temp;
    end
end
