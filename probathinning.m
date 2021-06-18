% obliczanie teoretycznych prawdopodobieńst i epirycznych wyników dla 
% algorytmu thinning
l=[0.2;0;0.6;0.7;0.25;0.9;0.35];
v=orth(hilb(7));
K=v*diag(l)*v';
L=K/(eye(7)-K);
pusty=[];
prawiepusty=[1,3];
sredni=[1,2,4,6];
prawiepelny=[1,2,3,5,6,7];
ile1=0;
ile2=0;
ile3=0;
ile4=0;
n=15000;
for i=1:n
    Y=thinning(K);
    if size(Y,2)==size(pusty)
        ile1=ile1+1;
    end
    if size(Y)==size(prawiepusty)
        if Y==prawiepusty 
            ile2=ile2+1;
        end
    end
    if size(Y)==size(sredni)
        if Y==sredni
            ile3=ile3+1;
        end
    end
    if size(Y)==size(prawiepelny)
        if Y==prawiepelny
            ile4=ile4+1;
        end
    end
end
d=det(L+eye(7));
pr1=1/d
disp(ile1/n)
pr2=det(L(prawiepusty,prawiepusty))/d
disp(ile2/n)
pr3=det(L(sredni,sredni))/d
disp(ile3/n)
pr4=det(L(prawiepelny,prawiepelny))/d
disp(ile4/n)
