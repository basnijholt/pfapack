%some tests on real matrices

A=rand(8,8);
A=A-A.';

%All of the following should give the same result
disp 'Compute the Pfaffian of a test matrix using different methods'
disp 'Pfaffian from LTL decomposition:' 
pfaffian_LTL(A)
disp 'Pfaffian from Householder tridiagonalization:'
pfaffian_householder(A)
disp 'Pfaffian from Hessenberg decomposition:'
pfaffian_hessenberg(A)
disp 'All of the results above should be the same'
disp '-------------------------------------------'

%compare to determinant
disp 'Compare square of Pfaffian to determinant'
disp 'Pfaffian squared: ', pfaffian_LTL(A)^2
disp 'Determinant: ', det(A)
disp 'Should again be identical'
disp '-------------------------'

%Test decompositions
disp 'Compute the various decompositions'
disp 'Error of A(P,P)-L*T*transpose(L):'
[T,L,P]=skew_LTL(A);
norm(A(P,P)-L*T*L.')
disp 'The error should be small (<1e-14)'
disp 'Error of A-Q*T*transpose(Q):'
[T,Q]=skew_tridiagonalize(A);
norm(A-Q*T*Q.')
disp 'The error should be small (<1e-14)'
disp '----------------------------------'