tspan = 40;
n = 2;
tsample = 0.01;
x0 = [7;10];

phi = eye(2);
Gamma = tsample*eye(2);

K_obs = 0.05*[-90 100;
        -70 -80];

l = -K_obs*[-4.55; -6.36];

X_obs = zeros(n,tspan/tsample);

X_obs(:,1) = x0;

for i =1:tspan/tsample - 1
    X_obs(:,i+1) = (phi + Gamma*K_obs)*X_obs(:,i) + Gamma*l;  
end

disp(eig(phi + Gamma*K_obs))
disp(norm(0.9250 + 0.0497i))
plot(X_obs(1,:),X_obs(2,:))
disp(X_obs(:,tspan/tsample))