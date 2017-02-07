function [mkernel,mkernel_sum] = get_mkernel(lambda,r,bound,multi,varargin)
%GET_MKERNEL This function returns the mkernel and the Renyi entropy for a vector of random
%values r 
%   lambda --> lambda parameter asociated with r for kernel processing
%   sigma_k --> sigma parameter asociated with r for kernel processing
%   r --> vector of random values
%   multi --> indicates dimensionality, for multi > 1, the mkernel will be
%   calculated by blocks of multi

    if isempty(varargin)
        P = 1:1:length(r);
    else
        P = varargin{1};
    end
    
    if multi < 2
        mkernel = zeros(bound(1,2)-bound(1,1)+1,bound(2,2)-bound(2,1)+1);
        temp = 0;
        for i = bound(1,1):bound(1,2)
            for j = bound(2,1):bound(2,2)
                mkernel(i-bound(1,1)+1,j-bound(2,1)+1) = exp(-lambda*((r(P(j))-r(P(i))).^2));
                temp = mkernel(i-bound(1,1)+1,j-bound(2,1)+1)+temp;
            end
        end
        mkernel_sum = temp;
        
    else
        
        
        N = length(r);
        s = zeros(N-multi+1,multi);
        for i = 1:N-multi+1
            for j = 1:multi
                s(i,j) = r(i+j-1);
            end
        end
            
        if (bound(1,2)-bound(1,1) < multi) && (bound(2,2)-bound(2,1) >= multi)
            mkernel = zeros(bound(1,2)-bound(1,1)+1,bound(2,2)-bound(2,1)+1-multi+1);
            temp = 0;
            for i = bound(1,1):bound(1,2)
                for j = bound(2,1):bound(2,2)-multi+1
                    mkernel(i-bound(1,1)+1,j-bound(2,1)+1) = exp(-lambda*(norm(s(P(j),:)-s(P(i),:))^2));
                    temp = mkernel(i-bound(1,1)+1,j-bound(2,1)+1)+temp;
                end
            end
            mkernel_sum = temp;
            
        elseif (bound(2,2)-bound(2,1) < multi) && (bound(1,2)-bound(1,1) >= multi)
            
            mkernel = zeros(bound(1,2)-bound(1,1)+1-multi+1,bound(2,2)-bound(2,1)+1);
            temp = 0;
            for i = bound(1,1):bound(1,2)-multi+1
                for j = bound(2,1):bound(2,2)
                    mkernel(i-bound(1,1)+1,j-bound(2,1)+1) = exp(-lambda*(norm(s(P(j),:)-s(P(i),:))^2));
                    temp = mkernel(i-bound(1,1)+1,j-bound(2,1)+1)+temp;
                end
            end
            mkernel_sum = temp;
            
        else
            
            mkernel = zeros(bound(1,2)-bound(1,1)+1-multi+1,bound(2,2)-bound(2,1)+1-multi+1);
            temp = 0;
            for i = bound(2,1):bound(2,2)-multi+1
                for j = bound(1,1):bound(1,2)-multi+1
                    mkernel(i-bound(1,1)+1,j-bound(2,1)+1)= exp(-lambda*(norm(s(P(j),:)-s(P(i),:))^2));
                    temp = mkernel(i-bound(1,1)+1,j-bound(2,1)+1)+temp;
                end
            end
            mkernel_sum = temp;
        end
    end
        
        
        
        
        
        
%         
%         mkernel = zeros(bound(1,2)-bound(1,1)+1,bound(2,2)-bound(2,1)+1);
%         temp = 0;
%         for i = bound(1,1):bound(1,2)
%             for j = bound(2,1):bound(2,2)
%                 mkernel(i-bound(1,1)+1,j-bound(2,1)+1) = exp(-lambda*(norm(s(P(j),:)-s(P(i),:))^2));
%                 temp = mkernel(i-bound(1,1)+1,j-bound(2,1)+1)+temp;
%             end
%         end
%         mkernel_sum = temp;


end

