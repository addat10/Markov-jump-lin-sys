function L=gen_rand_Laplacian(N,density)
A=rand(N)<density;
A=(A+A');
A(find(A>1))=1;
L=diag(A*ones(N,1))-A;
end