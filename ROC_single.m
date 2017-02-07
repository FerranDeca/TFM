n = 2;
N = 2000; %points for each block
M = 200; %number of blocks
multi = 1;
method = struct;
method.name = 'cholesky'; % 'kernel_matrix' for normal computation, 'cholesky' 
%                           for using cholesky decomposition
method.error = 0.1; %for cholesky method only


%% Generate and measure dependent variables

mu = [0 0];
rho = 0.3; %max 1/sqrt(2)
r = generate_gaussian_samples(N*M,mu,rho);
sign = round(rand(N*M,1));
sign( sign==0 )=-1;
r(:,1) = r(:,1).*sign;

sigma = [var(r(:,1)) var(r(:,2))];
figure;
plot(r(:,1),r(:,2),'o')
title('Generated dependent random values')

I_2_1 = zeros(M,1);

datestr(now)
disp('Començament');
% figure;
% try
    for i = 1:M
    I_2_1(i,1) = I2_joint_Parzen(r((i-1)*(N)+1:i*N,:),sigma,method,multi); 
    end
%     hold off;
% catch
%     disp('method.name is empty o not correct, or using cholesky would not gain anything, use kernel_matrix as the method name');
%     return
% end
datestr(now)
disp('Acaba mètode');


%% generate and measure independent variables

mu = [0 0];
rho = 0;
r2 = generate_gaussian_samples(N*M,mu,rho);

figure;
plot(r2(:,1),r2(:,2),'o')
title('Generated independent random values')

I_2_2 = zeros(M,1);
% figure;
for i = 1:M
   I_2_2(i,1) = I2_joint_Parzen(r2((i-1)*(N)+1:i*N,:),sigma,method,multi); 
end
% hold off;
datestr(now)
disp('Acaba mesures');


%% Draw ROC

roc = get_roc(10000,95,I_2_1,I_2_2);
figure; draw_roc(roc);
figure; stem(I_2_1); hold on; stem(I_2_2);
 
datestr(now)


