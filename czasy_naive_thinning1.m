% porównanie czasów działania dla algorytów naive i thinning w zależności
% od wymiaru przestrzeni
q=[100,500,1000,1500,2250,3000];
N=zeros(1,6);
T=zeros(1,6);
a=1;
for i=q
    l=rand(i,1);
    l=l/sum(l)*20;
    v=orth(hilb(i)+eye(i));
    K=v*diag(l)*v';
    f=@()naive(K);
    g=@()thinning(K);
    n=timeit(f);
    N(a)=n;
    t=timeit(g);
    T(a)=t;
    a=a+1;
end
plot(q,N,q,T,'Marker','*','MarkerSize',10,'LineWidth',2);
legend('Algorytm sekwencyjny','Algorytm przerzedzający','Location','northwest');
xlabel('Liczność przestrzeni');
ylabel('czas [s]');