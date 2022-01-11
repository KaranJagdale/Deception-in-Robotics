Kobs = 0.1*[-1  -4;2 -3];
lobs = [-3;-1];
x0 = [7;10];

T_int = 0.1; %values at each second 
T_span = 40; %Total time of trajectory generation
Xobs = x0;
for i = 1:T_span/T_int
    
    [t,X] = ode45(@(t,x)dyn(t,x,Kobs,lobs),[(i-1)*T_int, i*T_int],Xobs(:,size(Xobs,2)));
    Xobs = [Xobs X(size(X,1),:)'];
end
%[t,X] = ode45(@(t,x)dyn(t,x,Kobs,lobs),[0,1000],x0);
plot(Xobs(1,:),Xobs(2,:))
%plot(X(:,1),X(:,2));
disp(-Kobs\lobs)
