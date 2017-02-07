
datestr(now)
n = 2;
N_final = 2000;
methods = 4;
roc_points = 1000;
M = 300;
loops = 200;
multi = 1;
method = struct;
method.name = 'cholesky'; % 'kernel_matrix' for normal computation, 'cholesky' 
%                           for using cholesky decomposition
method.error = 0.1; %for cholesky method only

%% Generate dependent (r) and independent (r2) random values

mu = [0 0];
rho = 0.3;
r = generate_gaussian_samples(N_final*M,mu,rho);
sign = round(rand(N_final*M,1));
sign( sign==0 )=-1; 
r(:,1) = r(:,1).*sign;
figure;
plot(r(:,1),r(:,2),'o')
title('Generated dependent random values')
sigma1 = [var(r(:,1)) var(r(:,1))];

mu = [0 0];
rho = 0;
r2 = generate_gaussian_samples(N_final*M,mu,rho);
figure;
plot(r2(:,1),r2(:,2),'o')
title('Generated independent random values')
sigma2 = sigma1;

%% Loop until getting desired point

done = 0;
flag = zeros(methods,1);

results = zeros(loops,3,methods);
roc_curves = zeros(roc_points,methods*2,loops);

for N = N_final/loops:N_final/loops:N_final
    
    I_2_1 = zeros(M,1,methods);
    I_2_2 = zeros(M,1,methods);
%     flag = zeros(methods,1); % Flag if it is desired to stop until some
%     metric is accomplished
    
    for i = 1:length(I_2_1)
        if flag(1) == 0
            I_2_1(i,1,1) = I2_Parzen(r((i-1)*N+1:i*N,:),sigma1,method,multi);
            I_2_2(i,1,1) = I2_Parzen(r2((i-1)*N+1:i*N,:),sigma2,method,multi); 
        end
    end
    for i = 1:length(I_2_1)
        if flag(2) == 0
            I_2_1(i,1,2) = I2_joint_Parzen(r((i-1)*N+1:i*N,:),sigma1,method,multi);
            I_2_2(i,1,2) = I2_joint_Parzen(r2((i-1)*N+1:i*N,:),sigma2,method,multi);
        end
    end
    for i = 1:length(I_2_1)
        if flag(3) == 0
            I_2_1(i,1,3) = I2_CS(r((i-1)*N+1:i*N,:),sigma1,method,multi);
            I_2_2(i,1,3) = I2_CS(r2((i-1)*N+1:i*N,:),sigma2,method,multi);
        end
    end
    for i = 1:length(I_2_1)
        if flag(4) == 0
            I_2_1(i,1,4) = I2_HS(r((i-1)*N+1:i*N,:),sigma1,method,multi);
            I_2_2(i,1,4) = I2_HS(r2((i-1)*N+1:i*N,:),sigma2,method,multi);
        end
    end
    
    
    for i = 1:methods
%         if flag(i) == 0  % if a flag is wanted, uncomment this
            
            
            roc = get_roc(roc_points,95,I_2_1(:,:,i),I_2_2(:,:,i));
            roc_curves(:,2*i-1:2*i,N/(N_final/loops)) = roc;
            
            % 1- AUC
            
            auc = trapz(roc(:,1),roc(:,2));

            % Euclidean distance;
%             d = 1;
%             for t = 1:length(roc(:,1))
%                 dist = sqrt(((roc(t,1))^2)+((1- roc(t,2))^2));
%                 if dist < d
%                     d = dist;
%                 end
%             end
            
            
            % Deflection;
            
            deflection = abs(mean(I_2_2(:,1,i))-mean(I_2_1(:,1,i)))^2;
            deflection = deflection/var(I_2_2(:,1,i));

            
            % Save results
            
            results(N/(N_final/loops),1,i) = 1- auc;
            results(N/(N_final/loops),2,i) = deflection;
            results(N/(N_final/loops),3,i) = N;
            
%         end  % if a flag is wanted, decoment this
    end

    datestr(now)
    disp(['Acaba loop amb ' int2str(N) ' punts']);
    
end

%% Plot results

x = N_final/loops:N_final/loops:N_final;
figure;

for i = 1:methods
    plot(x,results(:,1,i));
    hold on
end
title('1 - AUC');
legend('First detector','Second detector','Third detector','Fourth detector');
xlabel('N')
hold off


figure;
for i = 1:methods
    plot(x,results(:,2,i));
    hold on
end
title('Deflection');
legend('First detector','Second detector','Third detector','Fourth detector');
xlabel('N')
hold off


