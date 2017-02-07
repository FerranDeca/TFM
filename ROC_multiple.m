datestr(now)

n = 2;
N = 1000; %points for each block
M = 1000; %number of blocks
multi = 1;
method = struct;
method.name = 'cholesky'; % 'kernel_matrix' for normal computation, 'cholesky' 
%                           for using cholesky decomposition
method.error = 0.1; %for cholesky method only
plot_type = 'roc'; %'roc' for normal display or 'log' for the logarithmic display of the roc


%% Generate and measure dependent variables

mu = [0 0];
rho = 0.3;
r = generate_gaussian_samples(N*M,mu,rho);
sign = round(rand(N*M,1));
sign( sign==0 )=-1; 
r(:,1) = r(:,1).*sign;

sigma = [var(r(:,1)) var(r(:,2))];
figure;
plot(r(:,1),r(:,2),'o')
title('Generated dependent random values')

I_2_1 = zeros(M,4);

datestr(now)
disp('Començament dependent');

try
    for i = 1:M
       I_2_1(i,1) = I2_Parzen(r((i-1)*(N)+1:i*(N),:),sigma,method,multi); 
    end

    datestr(now)
    disp('Acaba mètode Parzen');

    for i = 1:M
       I_2_1(i,2) = I2_joint_Parzen(r((i-1)*(N)+1:i*(N),:),sigma,method,multi); 
    end

    datestr(now)
    disp('Acaba mètode Joint Parzen');

    for i = 1:M
       I_2_1(i,3) = I2_CS(r((i-1)*(N)+1:i*(N),:),sigma,method,multi) ;
    end

    datestr(now)
    disp('Acaba mètode Cauchy Schwarz');

    for i = 1:M
       I_2_1(i,4) = I2_HS(r((i-1)*(N)+1:i*(N),:),sigma,method,multi); 
    end
catch
    disp('method.name is empty o not correct');
end
datestr(now)
disp('Acaba mètode Hilbert Spaces i últim');


%% generate and measure independent variables

mu = [0 0];
rho = 0;
r = generate_gaussian_samples(N*M,mu,rho);

figure;
plot(r(:,1),r(:,2),'o')
title('Generated independent random values')

I_2_2 = zeros(M,4);

datestr(now)
disp('Començament independent');

for i = 1:M
   I_2_2(i,1) = I2_Parzen(r((i-1)*(N)+1:i*(N),:),sigma,method,multi); 
end

datestr(now)
disp('Acaba mètode Parzen');

for i = 1:M
   I_2_2(i,2) = I2_joint_Parzen(r((i-1)*(N)+1:i*(N),:),sigma,method,multi); 
end

datestr(now)
disp('Acaba mètode Joint Parzen');

for i = 1:M
   I_2_2(i,3) = I2_CS(r((i-1)*(N)+1:i*(N),:),sigma,method,multi); 
end

datestr(now)
disp('Acaba mètode Cauchy_Schwarz');

for i = 1:M
   I_2_2(i,4) = I2_HS(r((i-1)*(N)+1:i*(N),:),sigma,method,multi); 
end

datestr(now)
disp('Acaba mètode Hilbert Spaces i últim');


%% Draw ROC

figure;
for j = 1:4
    roc = get_roc(10000,95,I_2_1(:,j),I_2_2(:,j));
    draw_roc(roc,j,'roc'); hold on;
end
legend('First detector','Second detector','Third detector','Fourth detector');



figure; stem(I_2_1(:,1)); hold on; stem(I_2_2(:,1)); title('Parzen'); 
figure; stem(I_2_1(:,2)); hold on; stem(I_2_2(:,2)); title('Joint Parzen'); 
figure; stem(I_2_1(:,3)); hold on; stem(I_2_2(:,3)); title('Cauchy Schwarz'); 
figure; stem(I_2_1(:,4)); hold on; stem(I_2_2(:,4)); title('Hilbert Spaces'); 

datestr(now)


