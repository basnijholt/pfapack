function [pf] = pfaffian_householder(A)
%pfaffian_householder: Compute the Pfaffian of a skew-symmetric matrix
%
% pf = pfaffian_householder(A) computes the Pfaffian of the skew-symmetric matrix A
% using the Householder algorithm

    assert(ndims(A)==2, 'argument must be a matrix')
    assert(size(A,1)==size(A,2), 'argument is not skew-symmetric')
    %make sure input is skew-symmetric
    assert(norm(A+A.')<1e-14*size(A,1), 'argument does not seem skew-symmetric')
    
    N=size(A,1);
    
    if( mod(N,2) == 1 )
        pf = 0.0;
        return;
    end
    
    pf = 1.0;
    for i = 1:2:N-2
        %Find a Householder vector to eliminate the i-th column
        [v, tau, alpha] = gallery('house', A(i+1:N,i));

        %Note that alpha is the (i+1,i) entry, i.e. one must
        %multiply the Pfaffian with -alpha. Then, one must also
        %take into account the determinant of the Householder
        %reflection which is -1 if tau not equal to 0, and +1
        %if tau is equal to 0 (then the transformation is the identity)
        if( tau ~= 0.0 )
            pf = pf * alpha;
        else
            pf = pf * -alpha;
        end
        
        if( tau ~= 0.0 )
            %Update the matrix block T(i+1:N,i+1:N)
            w = tau*A(i+1:N, i+1:N)*conj(v);
            A(i+1:N,i+1:N) = A(i+1:N,i+1:N) + v*w.'-w*v.';
        end
    end
    
    pf = pf * A(N-1,N);
end
