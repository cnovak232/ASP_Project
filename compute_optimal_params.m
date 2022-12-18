function best_params = compute_optimal_params(x,xn,ref_noise,order)

for p = 1:length(order)
    mu = linspace(.00001,.1,300);
    ser_lms = zeros(1,length(mu));
    ser_nlms = zeros(1,length(mu));
    for i = 1:length(mu)
        xc_lms = perform_lms(xn,ref_noise,mu(i),order(p));
        xc_nlms = perform_nlms(xn,ref_noise,mu(i),order(p));
    
        ser_lms(i) = compute_ser(x,xc_lms);
        ser_nlms(i) = compute_ser(x,xc_nlms);
    end
%     
%     figure;
%     plot(mu,mse_lms,mu,mse_nlms);
%     xlabel('Step Size (mu)');
%     ylabel('MSE');
%     legend('LMS','NLMS');
    
    [min__lms,loc_lms] = max(ser_lms);
    best_mu_lms(p) = mu(loc_lms);
    [min_mse_nlms,loc_nlms] = max(ser_nlms);
    best_mu_nlms(p) = mu(loc_nlms);
    
    lambda = linspace(.98,1,10);
    ser_rls = zeros(1,length(lambda));
    for i = 1:length(lambda)
        xc_rls = perform_rls(xn,ref_noise,lambda(i),order(p));
    
        ser_rls(i) = compute_ser(x,xc_rls);
    end
%     figure;
%     plot(lambda,mse_rls);
%     xlabel('Forgetting Factor (lambda)');
%     ylabel('MSE');
%     legend('RLS');
     [min_mse_rls,loc_rls] = max(ser_rls);
     best_lam_rls(p) = lambda(loc_rls);
    
    gamma = linspace(.5,1,50);
    ser_afa = zeros(1,length(gamma));
    for i = 1:length(gamma)
        xc_afa = perform_afa(xn,ref_noise,gamma(i),order(p));
    
        ser_afa(i) = compute_ser(x,xc_afa);
    end
%     figure;
%     plot(gamma,mse_afa);
%     xlabel('Gain Step Size (gamma)');
%     ylabel('MSE');
%     legend('AFA');
    [min_mse_afa,loc_afa] = max(ser_afa);
    best_gam_afa(p) = gamma(loc_afa);
end

best_params = struct;
best_params.mu_lms = best_mu_lms;
best_params.mu_nlms = best_mu_nlms;
best_params.lam_rls = best_lam_rls;
best_params.gam_afa = best_gam_afa;

% figure;
% plot(order,best_mu_lms(1:17),order,best_mu_nlms(1:17), );
% xlabel('Filter Order');
% ylabel('Mean Square Error');
% title('Best Mu vs Filter Order');
% legend('LMS','NLMS');