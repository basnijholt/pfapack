function [pf] = pfaffian_hesseneberg(A)
%pfaffian_hessenberg: Compute the Pfaffian of a real skew-symmetric matrix
%
% pf = pfaffian_hessenberg(A) computes the Pfaffian of the skew-symmetric matrix A
% using the Hessenberg decomposition

    assert(ndims(A)==2, 'argument must be a matrix')
    assert(size(A,1)==size(A,2), 'argument is not skew-symmetric')
    assert(isreal(A), 'argument must be a real matrix')
    %make sure input is skew-symmetric
    assert(norm(A+A.')<1e-14*size(A,1), 'argument does not seem skew-symmetric')
    
    N=size(A,1);
    
    if( mod(N,2) == 1 )
        pf = 0.0;
        return;
    end
    
    [Q,H] = hess(A);
    
    pf = 1.0;
    for i=1:2:N-1
        pf = pf * H(i,i+1);
    end
    pf = pf * det(Q);
    
end
