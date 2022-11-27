%Part 1 
%Initial vector and constant veclocity 
clear;
x_0 = [-2;-1];
v_star = [0.5;0.7];
%Time step and sigma value
dt = 0.1;
sig = 0.07;
%Observations b_k
b = zeros(2,50);
%True positions 
x_star = zeros(2,50);
ni = 1500;
%Calculate b_k and x_star at different time points
for k = 1:ni
    t_k = dt*k;
    x_star(:,k) = x_0 + v_star*t_k;
    b(:,k) = funbk(x_star(:,k),v_star) + sig*randn(2,1);
end
    
%Part 2 
%no of particles 
n = 500;
%Initialize the observations 
b_0 = funbk(x_0,v_star) + sig*randn(2,1);
%Calculate the angular positions 
theta_0 = b_0(1) + sig*randn(n,1);
r_0 = 1.5 + .05*randn(n,1);
%Calculate initial positions 
x01 = (r_0.*cos(theta_0))';
x02 = (r_0.*sin(theta_0))';
x_0 = [x01; x02];
%Initial Velocity guess
v0bar = [2;5];
gamma = 0.3;
v_0 = zeros(2,n);
for i = 1:n
   v_0(:,i) = v0bar + gamma*randn(2,1);
end
%Initialize the weights
w_0 = (1/n).*ones(n,1);


%Part 3
xhat = zeros(2,n,ni);
vhat = zeros(2,n,ni);
eta = zeros(1,n,ni);
x = zeros(2,n,ni);
v = zeros(2,n,ni);
w = zeros(1,n,ni);

for j = 1:ni
    %Part a
    if j==1
        xhat(:,:,1) = x_0 + dt.*v_0;
        vhat(:,:,1) = v_0;
    else
        xhat(:,:,j) = x(:,:,j-1) + dt.*v(:,:,j-1);
        vhat(:,:,j) = v(:,:,j-1);
    end
   
   %Part b
   %Compute fitness weights 
   for k = 1:n
    temp1 = funbk(xhat(:,k,j),vhat(:,k,j));
    temp2 = norm(b(:,j)-temp1,2).^2;
    eta(1,k,j)= ((1./sig.^2).*temp2);
   end  

   %Normalize fitness weights
   min_eta = min(eta(1,:,j));
   eta(1,:,j) = eta(1,:,j)-min_eta;
   eta(1,:,j) = exp(-eta(1,:,j));

   eta(1,:,j) = eta(1,:,j)./(mean(eta(1,:,j))*n);



   %Part c 
   idx = [];
   for i = 1:n
       xk = discrete([1:n],eta(1,:,j));
       idx = [idx;xk];
   end



   %Part d
   alpha = 0.01;
   for k = 1:n
       x(:,k,j) = xhat(:,idx(k),j) + alpha.*randn(2,1);
       v(:,k,j) = vhat(:,idx(k),j) + alpha.*randn(2,1);
   end

  %Part e 
  for k = 1:n
      temp3 = funbk(x(:,k,j),v(:,k,j));
      temp4 = norm(b(:,j)-temp3,2).^2;
      temp5 = funbk(xhat(:,idx(k),j),vhat(:,idx(k),j));
      temp6 = norm(b(:,j)-temp5,2).^2;
      w(1,k,j) = (1./sig.^2).*(temp4-temp6);
  end
 
   min_w = min(w(1,:,j));
   w(1,:,j) = exp(-w(1,:,j)+min_w);
   w(1,:,j) = w(1,:,j)./norm(w(1,:,j),1);
   
end

idx = zeros(ni,1);
mvals = zeros(ni,1);
dists = zeros(ni,1);
for i = 1:ni
    [m ind] = max(w(1,:,i));
    mvals(i) = m;
    idx(i) =ind;
    dists(i) = norm(funbk(x(:,idx(i),i),v(:,idx(i),i))-b(:,i),2);
   
end


 z01 = [];
 z02 = [];
for i = 1:n
   z0=  funbk(x(:,i,1450),v(:,i,1450));
  
   z01 = [z01;z0(1)];
   z02 = [z02;z0(2)];


end








