function s = wilsoncowan_RK2_graddescent(tau, b, W, k, s_init, C, chunk, dt, burn, ...
    xtarget, tolerance, gamma)
%% Use Gradient Descent to optimize external input (s) so that average 
% firing rate x (averaged across time and ROIs) is xtarget +- tolerance
% Strategy:
% Simulate a chunk of data with initial guess of s
% Average across time and space
% Set new s: s_(n+1) = s_n - gamma (<x>-xtarget)
% gamma is the learning rate
% Continue the procedure when error is within tolerance bounds, N times (5 maybe?)
% If the error doesn't go out of the tolerance bounds, stop
arguments
    tau                    % timescale
    b                      % bias term
    W                      % connectivity matrix
    k                      % sharpness of f-I curve
    s_init                 % input
    C                      % coupling strength
    chunk = 50             % time (s)
    dt = 10 / 1000         % time step (s)
    burn = 20              % burnin period (s)                
    xtarget = 0.1          % Targeted <x>
    tolerance = 0.02       % xtarget +- tolerance
    gamma = 0.1;           % Learning rate
end

N = 5; % I don't want to make this another adjustable parameter. Just set to 5 here. 
eps = 9999; % Starting error
s = s_init;
c = 1;
while true
    x = wilsoncowan_RK2(tau, b, W, k, s, C, chunk, dt, burn);
    meanx = mean(x(:));
    eps = meanx - xtarget;
    s = s - gamma * eps;
    % disp(eps)
    % disp(s)
    if abs(eps) < tolerance
        for i = 1:N
            x = wilsoncowan_RK2(tau, b, W, k, s, C, chunk, dt, burn);
            meanx = mean(x(:));
            eps = meanx - xtarget;
            if abs(eps) >= tolerance
                flag = false; % ;_;
                break
            else
                flag = true;  % :3
                continue
            end
        end
        if flag
            break
        end
    end
    % for debugging
    % subplot(1,2,1)
    % scatter(c, eps, 'k', 'filled'), hold on
    % subplot(1,2,2)
    % scatter(c, s, 'k', 'filled'), hold on
    c = c + 1;
    if c > 500
        error(sprintf("Doesn't converge, change stuff.\n Final s = %.3f, final eps = %.3f", s, eps))
    end
end
fprintf("\nFinal s = %.3f, final eps = %.3f", s, eps)
end
