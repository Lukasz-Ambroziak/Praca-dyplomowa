q=[100,500,1000,1500,2250,3000];
S=zeros(1,6);
T=zeros(1,6);
a=1;
for i=q
    l=rand(i,1);
    l=l/sum(l)*20;
    w=zeros(i,1);
    for j=1:i
    w(j)=l(j)/(1-l(j));
    end
    v=orth(hilb(i)+eye(i));
    K=v*diag(l)*v';
    f=@()spect(w,v);
    g=@()thinning(K);
    s=timeit(f);
    S(a)=s;
    t=timeit(g);
    T(a)=t;
    a=a+1;
end
plot(q,S,'g',q,T,'Marker','*','MarkerSize',10,'LineWidth',2);
legend('Algorytm spektralny','Algorytm przerzedzający','Location','northwest');
xlabel('Liczność przestrzeni');
ylabel('czas [s]');