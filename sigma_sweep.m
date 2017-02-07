datestr(now)

n = 2;
N = 500;
block = 500;
sigma = zeros(2,2);
roc_points = 1000;
sigma_steps = 0.5:0.5:2;
multi = 3;
method = struct;
method.name = 'kernel_matrix'; % 'kernel_matrix' for normal computation, 'cholesky' 
%                           for using cholesky decomposition
method.error = 0.1; %for cholesky method only

%% Generate and measure dependent variables

mu = [0 0];
rho = 0.3;
r = generate_gaussian_samples(N*block,mu,rho);
sign = round(rand(N*block,1));
sign( sign==0 )=-1; 
r(:,1) = r(:,1).*sign;

figure;
plot(r(:,1),r(:,2),'o')
title('Generated dependent random values')

sigma(1,:) = [var(r(:,1)) var(r(:,1))];


%% generate and measure independent variables

mu = [0 0];
rho = 0;
r2 = generate_gaussian_samples(N*block,mu,rho);

figure;
plot(r2(:,1),r2(:,2),'o')
title('Generated independent random values')

sigma(2,:) = sigma(1,:);

%% Get roc for every loop



for j = 1:length(sigma_steps)
    
    I_2_1 = zeros(block,1);
    I_2_2 = zeros(block,1);
    

    for i = 1:block
        I_2_1(i,1) = I2_CS(r((i-1)*(N)+1:i*(N),:),sigma(1,:)*sigma_steps(j),method,multi);
        I_2_2(i,1) = I2_CS(r2((i-1)*(N)+1:i*(N),:),sigma(2,:)*sigma_steps(j),method,multi); 
    end
    
    legend_info{j} = [num2str(sigma_steps(j)) ' \sigma_1'];
    roc = get_roc(roc_points,95,I_2_1(:,1),I_2_2(:,1));
    draw_roc(roc,j); hold on;

    
    datestr(now)
    disp(['Acaba ' int2str(j) ' loop de ' int2str(length(sigma_steps))]);
end

legend(legend_info);







