function [roc] = get_roc(N,percentile,I_2_1,I_2_2)
%GET_ROC returns the roc curve values in a Nx2 matrix, each column
%representing the y and x axis.
%   N --> Number of values to calculate (typical 10000)
%   percentile --> the values outside this percentile wil not take account
%   on roc measure (typical value 95)
%   I_2_1 and I_2_2 --> information measures for the dependent and
%   independent hypotesis

threshold = zeros(N,1);

roc = zeros(length(threshold),2);


tope1 =  prctile(I_2_1,percentile);
tope2 =  prctile(I_2_2,percentile);

if tope2> tope1
    for t = 1:length(threshold)
        step = tope2/length(threshold);
        threshold(t) = step*(t-1);
    end
else
    for t = 1:length(threshold)
        step = tope1/length(threshold);
        threshold(t) = step*(t-1);
    end  
end

for t = 1:length(threshold)
    sens = 0;
    fpositive = 0;
    for i = 1:length(I_2_2)
        if threshold(t) >= I_2_1(i)
            sens = sens + 1;
        end
        if threshold(t) >= I_2_2(i)
            fpositive = fpositive+1;
        end
    end
    roc(t,1) = sens/length(I_2_2);
    roc(t,2) = fpositive/length(I_2_2);
end

end

