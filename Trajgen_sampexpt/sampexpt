Km = 4; %max frobenious norm of matrix K
lm = 10; %max norm of l
n  = 2; %dim of K
n_samp = 50000;
samp_t = 0.1; %Sampling time of the discretised system 
xt = [-10 -4.55;
      0 -6.36];
nxt = 2;  %number of targets
phi = eye(2);
Gamma = samp_t*eye(2);
eps = 10;  %epsilon in the inequality below (28)
n_obs = 30; %number of observations
%To sample uniformy from above K and l we'll use rejection method

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
     if max(abs(eig(KfM))) < 1         
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
dec = 1;
I_f = {};  %indices set satifying all the conditions

for i = 1:nxt
    I_f{i} = [];
end

count = 0;

% for j = 1:nxt
%     for i = 1:size(I{j},2)
%         for k =2 : n_obs
%             
%             ilhs = norm(state((k-2),vectom(Kxt(:,I{j}(i)),n), phi, Gamma, lxt{j}(:,i), Xobs(:,k-1)) - Xobs(:,k)); %factor of 10 because the difference in the sample times
%            
%             if ilhs > eps
%                  dec = dec*0;
%             end
%         end
%         if dec ==1
%             I_f{j} = [I_f{j} i];
%         end
%     end
% end

for j = 1:nxt
    for i = 1:size(I{j},2)
        for k =1 : n_obs-1
            
            count = count +1;
            ilhs = norm(Xobs(:,k+1) - (phi + Gamma*vectom(Kxt(:,I{j}(i)),n))*Xobs(:,k) - Gamma*lxt{j}(:,i));

            if ilhs > (1 + norm(phi + Gamma*vectom(Kxt(:,I{j}(i)),n)))*eps
                 dec = dec*0;
            end
        end
        if dec ==1
            I_f{j} = [I_f{j} i];
        end
    end
end
%disp(count)
%disp(size(I{1}))

% disp(state(0,vectom(Kxt(:,I{1}(1)),n),  phi, Gamma, lxt{j}(:,1), Xobs(:,1)))
% disp(state(1,vectom(Kxt(:,I{1}(1)),n),  phi, Gamma, lxt{j}(:,1), Xobs(:,1)))
% disp(state(2,vectom(Kxt(:,I{1}(1)),n),  phi, Gamma, lxt{j}(:,1), Xobs(:,1)))
% disp(Xobs(:,2))
% disp(Xobs(:,3))
% disp(Xobs(:,4))


