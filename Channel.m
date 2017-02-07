n = 2;
N = 500; %points for each block
M = 2000; %number of blocks
multi = 4;
method = struct;
method.name = 'kernel_matrix'; % 'kernel_matrix' for normal computation, 'cholesky' 
%                           for using cholesky decomposition
method.error = 0.1; %for cholesky method only
channel = 'Hammerstein';

%% Generate data

mu = 0;
variance = 1;
x =  mvnrnd(mu,variance,N*M);

if strcmp(channel,'Memoryless')
    
    r_1 = mvnrnd(mu,variance,N*M);
    
    % non-linear transformation
    sign = round(rand(N*M,1));
    sign( sign==0 )=-1; 
    r_2 = r_1.*sign;
    
    %noise
    sigma_noise = 2;
    noise = mvnrnd(0, sigma_noise, length(r_2));
    r_4=r_2+noise;
    r = [r_1 r_4];

elseif strcmp(channel,'Hammerstein')
    
    h = [0.2294 0.4588 0.6882 0.4588 0.2294];
    r_1 = mvnrnd(mu,variance,N*M);
    
    % non-linear transformation
    sign = round(rand(N*M,1));
    sign( sign==0 )=-1; 
    r_2 = r_1.*sign;
 
    %filter
    r_3 = conv(r_2,h);
    
    %noise
    sigma_noise = 0.5;
    noise = mvnrnd(0, sigma_noise, length(r_3));
    r_4=r_3+noise;
    r = [r_1 r_4(1:N*M)];
    
elseif strcmp(channel,'Wienner')  
    
    h = [0.2294 0.4588 0.6882 0.4588 0.2294];
    r_1 = mvnrnd(mu,variance,N*M);
    
    %filter
    r_2 = conv(r_1,h);
    
    % non-linear transformation
    sign = round(rand(N*M,1));
    sign( sign==0 )=-1; 
    r_3 = r_2(1:N*M).*sign;
    
    %noise
    sigma_noise = 2;
    noise = mvnrnd(0, sigma_noise, length(r_3));
    r_4=r_3+noise;
    r = [r_1 r_4];

    
elseif strcmp(channel,'MF')
    
    h = [0.2294 0.4588 0.6882 0.4588 0.2294];
    r_1 = mvnrnd(mu,variance,N*M);
    
    % non-linear transformation
    sign = round(rand(N*M,1));
    sign( sign==0 )=-1; 
    r_2 = r_1.*sign;
 
    %filter
    r_3 = conv(r_2,h);

    %noise
    sigma_noise = 2;
    noise = mvnrnd(0, sigma_noise, length(r_3));
    r_4=r_3+noise;
    
    %matched filter
    r_5 = conv(r_4,h);
    r_4 = r_5;
    r = [r_1 r_4(1:N*M)];
    
    
end




%% estimate

I_2_1 = zeros(M,1);

datestr(now)
disp('Començament');
sigma = [var(r_1) var(r_1)];
% try
    for i = 1:M
    I_2_1(i,1) = I2_CS(r((i-1)*(N)+1:i*N,:),sigma,method,multi); 
    end
% catch
%     disp('method.name is empty o not correct, or using cholesky would not gain anything, use kernel_matrix as the method name');
%     return
% end
datestr(now)
disp('Acaba mètode');

I_2_2 = zeros(M,1);
r2 = [x r_4(1:N*M)];

for i = 1:M
   I_2_2(i,1) = I2_CS(r2((i-1)*(N)+1:i*N,:),sigma,method,multi); 
end

datestr(now)
disp('Acaba mesures');


%% Draw ROC

roc = get_roc(10000,95,I_2_1,I_2_2);
figure; draw_roc(roc);
% figure; stem(I_2_1); hold on; stem(I_2_2);
 
datestr(now)


