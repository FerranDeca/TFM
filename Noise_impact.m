n = 2;
N = 500; %points for each block
M = 2000; %number of blocks
multi = 1;
method = struct;
method.name = 'cholesky'; % 'kernel_matrix' for normal computation, 'cholesky' 
%                           for using cholesky decomposition
method.error = 0.1; %for cholesky method only

mu = 0;
variance = 1;
x =  mvnrnd(mu,variance,N*M);
r_1 = mvnrnd(mu,variance,N*M);
    % non-linear transformation
    sign = round(rand(N*M,1));
    sign( sign==0 )=-1; 
    r_2 = r_1.*sign;

figure;
for j = 1:5
    
    %noise
    sigma_noise = 1*j;
    noise = mvnrnd(0, sigma_noise, length(r_2));
    r_4=r_2+noise;
    r = [r_1 r_4];
    
    I_2_1 = zeros(M,1);

    sigma = [var(r_1) var(r_1)];

        for i = 1:M
        I_2_1(i,1) = I2_CS(r((i-1)*(N)+1:i*N,:),sigma,method,multi); 
        end


    I_2_2 = zeros(M,1);
    r2 = [x r_4(1:N*M)];

    for i = 1:M
       I_2_2(i,1) = I2_CS(r2((i-1)*(N)+1:i*N,:),sigma,method,multi); 
    end


    roc = get_roc(10000,95,I_2_1,I_2_2);
    draw_roc(roc,j);
    hold on;
    
end

legend('SNR = 0dB', 'SNR = -3dB','SNR = -6dB','SNR = -9dB','SNR = -12dB');

