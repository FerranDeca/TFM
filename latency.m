
N = 2000;
mu = 0;
sigma = 1;
sigma_noise = 0.1;
multi = 4;
real_delay = N/2;
block_size = N/4;
method = struct;
method.name = 'kernel_matrix';
method.error = 0.1;

%% Generate and measure dependent variables

r = mvnrnd(mu,sigma,N);

% figure;
% plot(r)
% title('Generated random values');

%% Process signal: 1 to 1 transform

r2 = r(real_delay:real_delay+block_size);
signe = round(rand(length(r2),1));
signe( signe==0 )=-1;
r2 = r2.*signe;

% figure;
% plot(r2)
% title('1 to 1 transform');

%% Filtering


h = [0.2294 0.4588 0.6882 0.4588 0.2294];
r3 = conv(r2,h);
% figure;
% plot(r3)
% title('Filtering');

%% Noise addition

noise = mvnrnd(0, sigma_noise, block_size+length(h));
r4 = noise + r3;
% figure;
% plot(r4)
% title('Noise addition');





%% Measure

m = 0;
sigma_kernel = [var(r) var(r)];
results = zeros(N-length(h)+1-block_size,1);
for i = 1:length(results)
    finestra = r(i:length(r4)+i-1);
    measure = I2_CS([finestra r4],sigma_kernel,method,multi);
    results(i) = measure;
    if measure > m
        m = measure;
        delay = i;
    end
end

display(delay);
display(real_delay);
n = 1:1:N-block_size-length(h)+1;
stem(n,results);






