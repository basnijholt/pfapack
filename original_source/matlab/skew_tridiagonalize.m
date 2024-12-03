function [T, Q] = skew_tridiagonalize(A)
%skew_tridiagonalize: Compute the tridiagonal form of A under unitary
%congurence (orthogonal similarity, if A is real)
%
% [T,Q] = skew_tridiagonalize(A) computes a skew-symmetric tridiagonal
% matrix T and a unitary Q such that A = Q*T*transpose(Q)

    assert(ndims(A)==2, 'argument must be a matrix')
    assert(size(A,1)==size(A,2), 'argument is not skew-symmetric')
    %make sure input is skew-symmetric
    assert(norm(A+A.')<1e-14*size(A,1), 'argument does not seem skew-symmetric')
    
    N=size(A,1);
    
    T=A;
    Q=eye(N);
    
    for i = 1:N-2
        %Find a Householder vector to eliminate the i-th column
        [v, tau, alpha] = gallery('house', T(i+1:N,i));
        
        T(i+1, i) = alpha;
        T(i, i+1) = -alpha;
        T(i+2:N, i) = 0;
        T(i, i+2:N) = 0;

        %Note: tau = 0 means the transformation is the identity
        if( tau ~= 0.0 )
            %Update the matrix block T(i+1:N,i+1:N)
            w = tau*T(i+1:N, i+1:N)*conj(v);
            T(i+1:N,i+1:N) = T(i+1:N,i+1:N) + v*w.'-w*v.';

            %Accumulate the individual Householder reflections
            %Accumulate them in the form P_1*P_2*..., which is
            %(..*P_2*P_1)^dagger
            y = tau*Q(:, i+1:N)* v;
            %Note: here really v', not v.' [v'=conj(v.')]
            Q(:, i+1:N)= Q(:, i+1:N) - y*v';
        end
    end
end

