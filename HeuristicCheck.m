%This code finds probability of orddering of heuristic volume being
%consistent with ordering of actual volume 
% i.e. V1<V2 and S1<S2 at the same time or V1>V2 and S1>S2 at the same time
group = 4;
n = 8;
degree = 4;
Tis = 1000;
count = 0;
ccount = 0; %complementary count
for k = 1:Tis
    ac = unifrnd(1,5,[2,8]);
    b = unifrnd(10,15,[1,8]);
    epsilon = 1;
    Root = [];
    for i = 1:n
        p = [ac(1,i) 0 -b(i) 0 (ac(2,i)-epsilon)];
        %disp(i)
        root =roots(p);
        root = sort(root, "descend");
        Root = [Root root];
    end    
    
    %Heuristic calculation
    Hroot = [];
    for i = 1:n/group
        hp = [0 0 0 0 0];
        for j = ((i-1)*group +1):i*group
            hp(1) = hp(1) + ac(1,j);
            hp(3) = hp(3) - b(j);
            hp(5) = hp(5) + ac(2,j) - epsilon;
        end 
        hroot = roots(hp);
        hroot = sort(hroot, "descend");
        Hroot = [Hroot hroot];
    end    
    sublevel = zeros(n/group, degree);
    for i  = 1:n/group
        for j = 1:degree
            if rem(j,2) == 1
                sublevel(i,j) = min(Root(j,((i-1)*group +1):i*group));
            end
            if rem(j,2) == 0
                sublevel(i,j) = max(Root(j,((i-1)*group +1):i*group));
            end
        end
    end
    %Here, volume is basically the interval between the roots as we are
    %dealing with polynomials. It is calculated as 2 times the difference
    %between largest and second largest root as we are dealing with quartic
    %polynomials and their roots are symmetric about 'y-axis'
    act1 = 2*(sublevel(1,1) - sublevel(1,2));  
    act2 = 2*(sublevel(2,1) - sublevel(2,2));
    heu1 = 2*(Hroot(1,1) - Hroot(2,1));
    heu2 = 2*(Hroot(1,2) - Hroot(2,2));
    
    %Comparing actual volume with heuristic volume
    if act1>act2 && heu1>heu2
        count = count +1;
    end    
    
    if (act1<act2) && (heu1<heu2)
        count = count +1;
    end    
    
    if act1<act2 && heu1>heu2
        ccount = ccount +1;
    end    
    
    if (act1>act2) && (heu1<heu2)
        ccount = ccount +1;
    end   
%     disp(['act1 : ', num2str(act1)])
%     disp(['act2 : ', num2str(act2)])
%     disp(['heu1 : ', num2str(heu1)])
%     disp(['heu2 : ', num2str(heu2)])
    
end
disp(count/Tis) %this is the probability that the ordering of heuristic and actual are consistent
disp(ccount/Tis) %Calculated this just to ensure that every case is covered
