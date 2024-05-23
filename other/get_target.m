function target = get_target(n,p)

if n==1000 && p==0.5, target = 5e-6;
elseif n==3000 && p==0.5, target = 5e-7;
elseif n==5000 && p==0.5, target = 5e-7;
elseif n==1000 && p==0.3, target = 5e-6;
elseif n==3000 && p==0.3, target = 5e-6;
elseif n==5000 && p==0.3, target = 5e-7;
elseif n==1000 && p==0.1, target = 0.2;
elseif n==3000 && p==0.1, target = 5e-3;
elseif n==5000 && p==0.1, target = 5e-5;
elseif n==1000 && p==0.05, target = 4;
elseif n==3000 && p==0.05, target = 0.4;
elseif n==5000 && p==0.05, target = 0.1;
else, target=5e-3; 
end
end

