function model_output = model_Kalman_Filter_matlab(params, free_choices, rewards, mdp)

    dbstop if error;
    G = mdp.G; % number of games
    T = mdp.T; % number of forced choices
    simmed_choices = nan(1,G); % initialize an empty array for simulated choices based on fit parameter values

    mu0 = 50; % initial reward value for a bandit is fixed at 50

    alpha_start = params.alpha_start; % starting learning rate parameter
    alpha_inf = params.alpha_inf; % asymptotic learning rate parameter
    info_bonuses = [params.info_bonus_h1 params.info_bonus_h6]; % information bonus parameters depending on horizon
    % decision noise parameters depending on horizon and information condition
    decision_noises = [params.dec_noise_h1_22 params.dec_noise_h6_22;  
                       params.dec_noise_h1_13 params.dec_noise_h6_13];
    % side bias parameters depending on horizon and information condition
    biases = [params.spatial_bis_h1_22 params.spatial_bis_h6_22;
              params.spatial_bis_h1_13 params.spatial_bis_h6_13];

    % calculate value for first learning rate      
    alpha0  = alpha_start / (1 - alpha_start) - alpha_inf^2 / (1 - alpha_inf);
    alpha_d = alpha_inf^2 / (1 - alpha_inf); 

    % initialize data structures
    action_probs = nan(1,G);
    pred_errors = nan(T+1,G);
    pred_errors_alpha = nan(T+1,G);
    exp_vals = nan(T+1,G);
    alpha = nan(T+1,G);
    
    for g=1:G  % loop over games
        % values
        mu1 = [mu0 nan nan nan nan];
        mu2 = [mu0 nan nan nan nan];

        % learning rates 
        alpha1 = [alpha0 nan nan nan nan]; 
        alpha2 = [alpha0 nan nan nan nan]; 

        % extract information bonus depending on horizon
        A = info_bonuses(mdp.horizon_sequence(g));
        % extract decision noise depending on horizon/information condition
        sigma_g = decision_noises(mdp.info_sequence(g), mdp.horizon_sequence(g));
        % extract side bias depending on horizon/information condition
        bias = biases(mdp.info_sequence(g), mdp.horizon_sequence(g));

        for t=1:T  % loop over forced-choice trials
            % if the first bandit was chosen
            if (mdp.forced_choices(t,g) == 1) 
                % update learning rate for each bandit
                alpha1(t+1) = 1/( 1/(alpha1(t) + alpha_d) + 1 );
                alpha2(t+1) = 1/( 1/(alpha2(t) + alpha_d) );
                exp_vals(t,g) = mu1(t);
                pred_errors(t,g) = (rewards(t,g) - exp_vals(t,g));
                alpha(t,g) = alpha1(t+1);
                % use learning rate and prediction error to updated
                % expected reward for the chosen bandit
                pred_errors_alpha(t,g) = alpha1(t+1) * pred_errors(t,g); % confirm that alpha here should be t+1
                mu1(t+1) = mu1(t) + pred_errors_alpha(t,g);
                mu2(t+1) = mu2(t); 
            else
                % update learning rate for each bandit
                alpha1(t+1) = 1/( 1/(alpha1(t) + alpha_d) ); 
                alpha2(t+1) = 1/( 1/(alpha2(t) + alpha_d) + 1 );
                exp_vals(t,g) = mu2(t);
                mu1(t+1) = mu1(t);
                pred_errors(t,g) = (rewards(t,g) - exp_vals(t,g));
                alpha(t,g) = alpha2(t+1);
                % use learning rate and prediction error to updated
                % expected reward for the chosen bandit
                pred_errors_alpha(t,g) = alpha2(t+1) * pred_errors(t,g);
                mu1(t+1) = mu1(t);
                mu2(t+1) = mu2(t) + pred_errors_alpha(t,g);
            end
        end
        % Aggregate the difference in utility between the right and left bandit based
        % on expected reward value, information bonus, and side bias
        dQ = mu2(T+1) - mu1(T+1) + A * mdp.right_info(g) + bias;

        % probability of choosing the right bandit
        p = 1 / (1 + exp(-dQ/(sigma_g)));

        action_probabilities = free_choices(g)*p + (1-free_choices(g))*(1-p);
        action_probs(g) = action_probabilities;

        % sample from the distribution for choosing the right bandit to
        % simulate behavior under these parameters
        u = rand(1,1);
        choice = u <= p;

        % store simulated behavior
        simmed_choices(g) = choice;


        
    end
    
    model_output.simmed_choices = simmed_choices;