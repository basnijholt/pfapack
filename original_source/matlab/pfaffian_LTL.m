function [ pf ] = pfaffian_LTL( A )
%pfaffian_LTL: Compute the Pfaffian of a skew-symmetric matrix
%
% pf = pfaffian_LTL(A) computes the Pfaffian of the skew-symmetric matrix A
% using the Parlett-Reid algorithm

    assert(ndims(A)==2, 'argument must be a matrix')
    assert(size(A,1)==size(A,2), 'argument is not skew-symmetric')
    %make sure input is skew-symmetric
    assert(norm(A+A.')<1e-14*size(A,1), 'argument does not seem skew-symmetric')

    N = size(A,1);
    
    if( mod(N,2)==1 )
        pf = 0.0;
        return
    end
    
    pf = 1.0;
    for k = 1:2:N-2
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
            
            %Every interchange gives a "-" for the determinant of the
            %permutation
            pf = -pf;
        end
        
        pf = pf * A(k,k+1);
        
        %Now form the Gauss vector
        if( A(k+1,k) ~= 0.0 )
            tau = A(k+2:N,k)/A(k+1,k);

            %Update the matrix block A(k+2:,k+2)
            A(k+2:N,k+2:N) = A(k+2:N,k+2:N) + tau*A(k+2:N,k+1).' - A(k+2:N,k+1)*tau.';
        end
    end

    pf = pf * A(N-1,N);
end

