% obliczanie teoretycznych prawdopodobieńst i epirycznych wyników dla
% algorytmu thinning  
l=orth(magic(100));
l=abs(l(:,1)*9);
v=orth(hilb(100)+eye(100));
K=v*diag(l)*v';
maly=[16,28,39,47,56,63,78,81];
prawiepusty=[3,7,15,19,23,29,34,39,45,58,67,73,88,91,95,99];
sredni=(1:2:100);
prawiepelny=setdiff((1:100),[12,34,51,58,72,80,88]);
ile1=0;
ile2=0;
ile3=0;
ile4=0;
n=15000;
for i=1:n
    Y=thinning(K);
    if ismember(maly,Y)==ones(1,8)
        ile1=ile1+1;
    end
    if ismember(prawiepusty,Y)==ones(1,16)
        ile2=ile2+1;
    end
    if ismember(sredni,Y)==ones(1,50)
        ile3=ile3+1;
    end
    if ismember(prawiepelny,Y)==ones(1,93)
        ile4=ile4+1;
    end
end
pr1=det(K(maly,maly))
disp(ile1/n)
pr2=det(K(prawiepusty,prawiepusty))
disp(ile2/n)
pr3=det(K(sredni,sredni))
disp(ile3/n)
pr4=det(K(prawiepelny,prawiepelny))
disp(ile4/n)