function [ T, L, P ] = skew_LTL( A )
%skew_LTL: Compute the L*T*transpose(L) decomposition of a skew-symmetric
%matrix
%
% [T, L, P] = skew_LTL(A) computes a skew-symmetric tridiagonal matrix T, a
% lower unit triangular matrix L and a permutation vector P such that
% A(P,P) = L*T*transpose(L). If the permutation matrix itself is desired,
% it can be computed as Pmat=eye(size(A,1)); Pmat=Pmat(P,:)

    assert(ndims(A)==2, 'argument must be a matrix')
    assert(size(A,1)==size(A,2), 'argument is not skew-symmetric')
    %make sure input is skew-symmetric
    assert(norm(A+A.')<1e-14*size(A,1), 'argument does not seem skew-symmetric')
    
    N = size(A,1);
    
    L = eye(N);
    Pv = 1:N;
    
    for k = 1:N-2
        %First, find the largest entry in A(k+1:N,k) and
        %permute it to A(k+1,k)
        [c, kp] = max(A(k+1:N, k));
        kp = kp + k;

        %Check if we need to pivot
        if( kp ~= k+1 )
            %interchange rows k+1 and kp
            temp = A(k+1,k:N);
            A(k+1,k:N) = A(kp,k:N);
            A(kp,k:N) = temp;

            %Then interchange columns k+1 and kp
            temp = A(k:N,k+1);
            A(k:N,k+1) = A(k:N,kp);
            A(k:N,kp) = temp;

            %permute L accordingly
            temp = L(k+1,1:k);
            L(k+1,1:k) = L(kp,1:k);
            L(kp,1:k) = temp;
            
            %accumulate the permutation vector
            temp = Pv(k+1);
            Pv(k+1) = Pv(kp);
            Pv(kp) = temp;
        end
        
        %Now form the Gauss vector
        if( A(k+1,k) ~= 0.0 )
            tau = A(k+2:N,k)/A(k+1,k);
     

            %clear eliminated row and column
            A(k+2:N,k) = 0.0;
            A(k,k+2:N) = 0.0;

            %Update the matrix block A(k+2:,k+2)
            A(k+2:N,k+2:N) = A(k+2:N,k+2:N) + tau*A(k+2:N,k+1).' - A(k+2:N,k+1)*tau.';

            L(k+2:N,k+1) = tau;
        end
    end

    T=A;
    P=Pv;
end

