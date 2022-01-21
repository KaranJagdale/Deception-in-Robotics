%rng('default')
%rng(1)

Km = 10; %max frobenious norm of matrix K
lm = 50; %max norm of l
n  = 2; %dim of K
n_samp = 1000000;
samp_t = 0.01; %Sampling time of the discretised system 
tspan = 40;
xt = [-2 -4.55;
      -8 -6.36];
nxt = 2;  %number of targets
phi = eye(2);
Gamma = samp_t*eye(2);
eps = 0.5;  %epsilon in the inequality below (28)
n_obs = 15; %number of observations
%To sample uniformly from above K and l we'll use rejection method

Kr = -Km + 2*Km*rand(n^2,n_samp);
lr = -lm + 2*lm*rand(n,n_samp);
Kf = [];
%lf = [];
for i=1:n_samp
    if norm(Kr(:,i)) < Km
        Kf = [Kf Kr(:,i)];
    end
end
Kxt = [];
lxt = {};

for i = 1:nxt   
    lxt{i} = [];
end

for i = 1: size(Kf,2)
     KfM = vectom(Kf(:,i),n);
     if max(abs(eig(phi+ Gamma*KfM))) < 1         
         Kxt = [Kxt Kf(:,i)];    %Ks satisfying the eigen condition
                %lxt = [lxt KfM*xt(:,j)]; %K and l converging to xt
         for j = 1 : nxt
            lxt{j} = [lxt{j} -KfM*xt(:,j)]; %K and l converging to xt
         end
     end    
end

%Checking the norm condition on l and store indices satifying the condition

I = {};   %indices set satisfying the l norm condition
for i = 1:nxt
    I{i} = [];
end

for j = 1:nxt
    for i= 1:size(lxt{j},2)
        if norm(lxt{j}(:,i)) < lm
            I{j} = [I{j} i];
        end
    end
end
dec = true;
I_f = {};  %indices set satifying all the conditions

for i = 1:nxt
    I_f{i} = [];
end

count = 0;

for j = 1:nxt
    for i = 1:size(I{j},2)
        for k =2 : n_obs
            
            ilhs = norm(state((k-2),vectom(Kxt(:,I{j}(i)),n), phi, Gamma, lxt{j}(:,I{j}(i)), X_obs(:,1)) - X_obs(:,k)); %factor of 10 because the difference in the sample times
           
            if ilhs > eps
                 dec = false;
                 continue
            end
        end
        if dec
            I_f{j} = [I_f{j} I{j}(i)];
        end
        dec = true;
    end
end
cond = 1;
% for j = 1:nxt
%     for i = 1:size(I{j},2)
%         for k =1 : n_obs-1
%             
%             count = count +1;
%             ilhs = norm(X_obs(:,k+1) - (phi + Gamma*vectom(Kxt(:,I{j}(i)),n))*X_obs(:,k) - Gamma*lxt{j}(:,I{j}(i)));
% 
%             if ilhs > (1 + norm(phi + Gamma*vectom(Kxt(:,I{j}(i)),n),'fro'))*eps
%                  dec = false; 
%                  continue
%             end            
%         end
%         if dec
%             I_f{j} = [I_f{j} I{j}(i)];
%         end
%         if dec && cond == 1
%             disp('Hey')
%             disp(norm((phi + Gamma*vectom(Kxt(:,I{j}(i)),n)),'fro')*eps + eps)
%             for z = 1:n_obs-1
%                 disp(norm(X_obs(:,z+1) - (phi + Gamma*vectom(Kxt(:,I{j}(i)),2))*X_obs(:,z) - Gamma*lxt{j}(:,I{j}(i))))
%             end 
%             cond =0;
%         end
%         dec = true;
%     end
% end

Kcl = [];
thres = 2;
for i = 1:tspan/samp_t
    if(norm(K_obs(1,1) - Kr(1,i))<thres) && (norm(K_obs(1,2) - Kr(2,i))<thres) && (norm(K_obs(2,1) - Kr(3,i))<thres) && (norm(K_obs(2,2) - Kr(4,i))<thres)
        Kcl = [Kcl Kr(:,i)];
    end
end

Xp = zeros(n,tspan/samp_t);
Kp = vectom(Kxt(:,I{1}(3)),n);
lp = lxt{1}(:,I{1}(3));
disp(lp)
disp(-Kp*xt(:,1))
e = eig(phi + Gamma*Kp);
disp(e)
disp(max(abs(eig(phi+Gamma*Kp))))
Xp(:,1) = X_obs(:,1);

Xp1 = zeros(n,tspan/samp_t);
Kp1 = vectom(Kxt(:,I{2}(2)),n);
lp1 = lxt{2}(:,I{2}(2));
disp(lp)
disp(-Kp*xt(:,2))
e = eig(phi + Gamma*Kp);
disp(e)
disp(max(abs(eig(phi+Gamma*Kp))))
Xp1(:,1) = X_obs(:,1);

for i = 1:(tspan/samp_t)-1
    Xp(:,i+1) = (phi + Gamma*Kp)*Xp(:,i) + Gamma*lp;
end

for i = 1:(tspan/samp_t)-1
    Xp1(:,i+1) = (phi + Gamma*Kp1)*Xp1(:,i) + Gamma*lp1;
end

Kg = vectom(Kxt(:,I_f{1}(1)),n);
lg = -Kg*xt(:,1);
Xg = zeros(n,tspan/samp_t);
Xg(:,1) = X_obs(:,1);
for i=1:tspan/samp_t -1
    Xg(:,i+1) = (phi + Gamma*Kg)*Xg(:,i) + Gamma*lg;
end

Kc = vectom(Kcl(:,1),n);
lc = -Kc*xt(:,2);
Xc = zeros(n,tspan/samp_t);
Xc(:,1) = X_obs(:,1);
for i=1:tspan/samp_t -1
    Xc(:,i+1) = (phi + Gamma*Kc)*Xc(:,i) + Gamma*lc;
end


figure(1)
plot(X_obs(1,:),X_obs(2,:))
hold on
plot(Xp(1,:),Xp(2,:))
plot(Xp1(1,:),Xp1(2,:))
plot(Xg(1,:),Xg(2,:))
plot(Xc(1,:),Xc(2,:))
legend('obs','test1', 'test2', 'cons', 'Xc')
%disp(count)
%disp(size(I{1}))
% figure(2)
% scatter3(Kr(4,:),Kr(2,:),Kr(3,:), '.')
% figure(3)
% scatter(Kr(1,:),Kr(2,:), '.')

Krfil = [];
for i = 1:tspan/samp_t
    if (Kr(3,i)) < 2
        Krfil = [Krfil Kr(:,i)];
    end
end

figure(4)
scatter(Krfil(1,:),Krfil(2,:), '.')

Kcl = [];
thres = 2;
for i = 1:tspan/samp_t
    if(norm(K_obs(1,1) - Kr(1,i))<thres) && (norm(K_obs(1,2) - Kr(2,i))<thres) && (norm(K_obs(2,1) - Kr(3,i))<thres) && (norm(K_obs(2,2) - Kr(4,i))<thres)
        Kcl = [Kcl Kr(:,i)];
    end
end


% disp(state(0,vectom(Kxt(:,I{1}(1)),n),  phi, Gamma, lxt{j}(:,1), Xobs(:,1)))
% disp(state(1,vectom(Kxt(:,I{1}(1)),n),  phi, Gamma, lxt{j}(:,1), Xobs(:,1)))
% disp(state(2,vectom(Kxt(:,I{1}(1)),n),  phi, Gamma, lxt{j}(:,1), Xobs(:,1)))
% disp(Xobs(:,2))
% disp(Xobs(:,3))
% disp(Xobs(:,4))


