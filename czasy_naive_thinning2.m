N=zeros(1,11);
T=zeros(1,11);
q=zeros(1,11);
v=orth(hilb(1000)+eye(1000));
for i=1:11
    if i>1
    l=rand(1000,1)*(i-1)/10;
    else
    l=rand(1000,1)/40;
    end
    q(i)=max(l); 
    K=v*diag(l)*v';
    f=@()naive(K);
    g=@()thinning(K);
    n=timeit(f);
    N(i)=n;
    t=timeit(g);
    T(i)=t;
end
plot(q,N,q,T,'Marker','*','MarkerSize',10,'LineWidth',2);
legend('Algorytm sekwencyjny','Algorytm przerzedzający','Location','northwest');
xlabel('maksymalna wartość własna');
ylabel('czas [s]');