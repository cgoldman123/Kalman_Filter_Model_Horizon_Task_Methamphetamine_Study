function fits = fit_model(formatted_file)
    rng(23);
    sub = load_data(formatted_file);

    disp(sub);
    
    %% ========================================================================
    %% HIERARCHICAL MODEL FIT ON FIRST FREE CHOICE %%
    %% HIERARCHICAL MODEL FIT ON FIRST FREE CHOICE %%
    %% HIERARCHICAL MODEL FIT ON FIRST FREE CHOICE %%
    %% ========================================================================

    %% prep data structure 
    L = unique(sub(1).gameLength); % length of games
    NS = length(sub);   % number of subjects
    T = 4;              % number of forced choices
    
    NUM_GAMES = max(vertcat(sub.game), [], 'all');

    a  = zeros(NS, NUM_GAMES, T); % forced choices (1=left bandit, 2=right bandit)
    c5 = nan(NS,   NUM_GAMES); % free choices (0=left bandit, 1=right bandit)
    r  = zeros(NS, NUM_GAMES, T); % rewards
    UC = nan(NS,   NUM_GAMES); % uncertainty condition (1=equal information, 2=unequal information)
    GL = nan(NS,   NUM_GAMES); % game length (1=Horizon of 1, 2=Horizon of 6)

    for sn = 1:length(sub)

        % choices on forced trials
        dum = sub(sn).a(:,1:4);
        a(sn,1:size(dum,1),:) = dum;

        % choices on free trial
        % note a slight hacky feel here - a is 1 or 2, c5 is 0 or 1.
        dum = sub(sn).a(:,5) == 2;
        L(sn) = length(dum);
        c5(sn,1:size(dum,1)) = dum;

        % rewards
        dum = sub(sn).r(:,1:4);
        r(sn,1:size(dum,1),:) = dum;

        % game length
        dum = sub(sn).gameLength;
        GL(sn,1:size(dum,1)) = dum;

        G(sn) = length(dum);

        % uncertainty condition 
        dum = abs(sub(sn).uc - 2) + 1;
        UC(sn, 1:size(dum,1)) = dum;

        % difference in information
        dum = sub(sn).uc - 2;
        dI(sn, 1:size(dum,1)) = -dum;
            

    end
    
    GL(GL==5) = 1;
    GL(GL==10) = 2;

    C1 = (GL-1)*2+UC;
    nC1 = 4;
    nC2 = 1;

    % meaning of condition 1 (SMT FIXED)
    % gl uc c1
    %  1  1  1 - horizon 1, [2 2]
    %  1  2  2 - horizon 1, [1 3]
    %  2  1  3 - horizon 6, [2 2]
    %  2  2  4 - horizon 6, [1 3]



    datastruct = struct(...
        'C1', C1, 'nC1', nC1, ...
        'NS', NS, 'G',  G,  'T',   T, ...
        'dI', dI, 'a',  a,  'c5',  c5, 'r', r);
    
    %% run hierarchical model fits! 
    nchains = 4;
    nburnin = 500;
    nsamples = 1000; 
    thin = 1;
    % MCMC parameters for JAGS


    % Initialize values all latent variables in all chains
    clear S init0
    for i=1:nchains

        S.a0(1:nC2) = 1;
        S.b0(1:nC2) = 1;
        S.a_inf(1:nC2) = 1;
        S.b_inf(1:nC2) = 1;
        S.AA(1:NS,1:nC1,1:nC2) = 0;
        S.BB(1:NS,1:nC1,1:nC2) = 100;
        init0(i) = S;
    end

    % Use JAGS to Sample
    tic
    doparallel = 1;

    fprintf( 'Running JAGS\n' );
    [samples, stats ] = matjags( ...
        datastruct, ...
        fullfile('./model_Kalman_Filter'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', thin, ...
        'monitorparams', ...
        {'a0' 'b0' 'alpha_start' 'alpha0' 'alpha_d' ...
        'a_inf' 'b_inf' 'alpha_inf' ...
        'mu0_mean' 'mu0_sigma' 'mu0'...
        'AA_mean' 'AA_sigma' 'AA'...
        'SB_mean' 'SB_sigma' 'SB' ...
        'BB_mean' 'BB' ...
        }, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0  );
    toc

    %% throw out first N-1 samples
    N = 1;

    stats.mean.SB = squeeze(mean(mean(samples.SB(:,N:end,:,:,:),2),1));
    stats.mean.BB = squeeze(mean(mean(samples.BB(:,N:end,:,:,:),2),1));
    stats.mean.AA = squeeze(mean(mean(samples.AA(:,N:end,:,:,:),2),1));

    %stats.mean.mu0 = squeeze(mean(mean(samples.mu0(:,N:end,:,:),2),1));
    stats.mean.alpha_start = squeeze(mean(mean(samples.alpha_start(:,N:end,:,:),2),1));
    stats.mean.alpha_inf = squeeze(mean(mean(samples.alpha_inf(:,N:end,:,:),2),1));
    stats.mean.alpha0= squeeze(mean(mean(samples.alpha0(:,N:end,:,:),2),1));
    stats.mean.alpha_d = squeeze(mean(mean(samples.alpha_d(:,N:end,:,:),2),1));
    
    
    
    %% Organize fits
    fits = struct();
    for si = 1:length({sub.subjectID})
        fits(si).id = {sub(si).subjectID};
      
        fits(si).info_bonus_h1 = stats.mean.AA(si, 2);
        fits(si).info_bonus_h6 = stats.mean.AA(si, 4);

        fits(si).dec_noise_h1_22 = stats.mean.BB(si, 1);
        fits(si).dec_noise_h1_13 = stats.mean.BB(si, 2);
        fits(si).dec_noise_h6_22 = stats.mean.BB(si, 3);
        fits(si).dec_noise_h6_13 = stats.mean.BB(si, 4);

        fits(si).spatial_bis_h1_22 = stats.mean.SB(si, 1);
        fits(si).spatial_bis_h1_13 = stats.mean.SB(si, 2);
        fits(si).spatial_bis_h6_22 = stats.mean.SB(si, 3);
        fits(si).spatial_bis_h6_13 = stats.mean.SB(si, 4);

        fits(si).alpha_start = stats.mean.alpha_start(si);
        fits(si).alpha_inf = stats.mean.alpha_inf(si);  
        
        % get variance
        % Extract samples after burn-in (N:end) and reshape to a vector
        % Compute variance of the samples

        person_samples = samples.AA(:, N:end, si, 2, 1);
        person_samples = person_samples(:);
        fits(si).info_bonus_h1_var = var(person_samples);
        
        person_samples = samples.AA(:, N:end, si, 4, 1);
        person_samples = person_samples(:);
        fits(si).info_bonus_h6_var = var(person_samples);
        
        person_samples = samples.BB(:, N:end, si, 1, 1);
        person_samples = person_samples(:);
        fits(si).dec_noise_h1_22_var = var(person_samples);
        
        person_samples = samples.BB(:, N:end, si, 2, 1);
        person_samples = person_samples(:);
        fits(si).dec_noise_h1_13_var = var(person_samples);      
        
        person_samples = samples.BB(:, N:end, si, 3, 1);
        person_samples = person_samples(:);
        fits(si).dec_noise_h6_22_var = var(person_samples);
        
        person_samples = samples.BB(:, N:end, si, 4, 1);
        person_samples = person_samples(:);
        fits(si).dec_noise_h6_13_var = var(person_samples);      
       
        person_samples = samples.alpha_start(:, N:end, si);
        person_samples = person_samples(:);
        fits(si).alpha_start_var = var(person_samples);
        
        person_samples = samples.alpha_inf(:, N:end, si);
        person_samples = person_samples(:);
        fits(si).alpha_inf_var = var(person_samples);
        
        
        params = fits(si);
        free_choices = c5(si,:);
        rewards = squeeze(r(si,:,:))';

        mdp.info_sequence = UC(si,:);
        mdp.horizon_sequence = GL(si,:);
        mdp.forced_choices = squeeze(a(1,:,:))';
        mdp.right_info = squeeze(dI(1,:,:));
        mdp.T = T; % num forced choices
        mdp.G = 80; % game length

        model_output(si).results = model_Kalman_Filter_matlab(params,free_choices, rewards, mdp);  

    end
    % Extract simulated free choices based on actual parameter fits
    for si = 1:length({sub.subjectID})
        simmed_c5(si, :) = model_output(si).results.simmed_choices;
    end
    
    %% Use JAGS to fit the simulated behavior


    datastruct = struct(...
    'C1', C1, 'nC1', nC1, ...
    'NS', NS, 'G',  G,  'T',   T, ...
    'dI', dI, 'a',  a,  'c5',  c5, 'r', r);


    datastruct.c5 = simmed_c5;


    fprintf( 'Running JAGS\n' );
        [simmed_samples, simmed_stats ] = matjags( ...
            datastruct, ...
            fullfile('./model_Kalman_Filter'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', thin, ...
            'monitorparams', ...
            {'a0' 'b0' 'alpha_start' 'alpha0' 'alpha_d' ...
            'a_inf' 'b_inf' 'alpha_inf' ...
            'mu0_mean' 'mu0_sigma' 'mu0'...
            'AA_mean' 'AA_sigma' 'AA'...
            'SB_mean' 'SB_sigma' 'SB' ...
            'BB_mean' 'BB' ...
            }, ...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0  );
        toc

    %% throw out first N-1 samples
    N = 1;

    simmed_stats.mean.SB = squeeze(mean(mean(simmed_samples.SB(:,N:end,:,:,:),2),1));
    simmed_stats.mean.BB = squeeze(mean(mean(simmed_samples.BB(:,N:end,:,:,:),2),1));
    simmed_stats.mean.AA = squeeze(mean(mean(simmed_samples.AA(:,N:end,:,:,:),2),1));

    %stats.mean.mu0 = squeeze(mean(mean(samples.mu0(:,N:end,:,:),2),1));
    simmed_stats.mean.alpha_start = squeeze(mean(mean(simmed_samples.alpha_start(:,N:end,:,:),2),1));
    simmed_stats.mean.alpha_inf = squeeze(mean(mean(simmed_samples.alpha_inf(:,N:end,:,:),2),1));
    simmed_stats.mean.alpha0= squeeze(mean(mean(simmed_samples.alpha0(:,N:end,:,:),2),1));
    simmed_stats.mean.alpha_d = squeeze(mean(mean(simmed_samples.alpha_d(:,N:end,:,:),2),1));

    %% Organize fits

     for si = 1:length({sub.subjectID})
        fits(si).simmed_id = {sub(si).subjectID};
    
        fits(si).simmed_info_bonus_h1 = simmed_stats.mean.AA(si, 2);
        fits(si).simmed_info_bonus_h6 = simmed_stats.mean.AA(si, 4);
    
        fits(si).simmed_dec_noise_h1_22 = simmed_stats.mean.BB(si, 1);
        fits(si).simmed_dec_noise_h1_13 = simmed_stats.mean.BB(si, 2);
        fits(si).simmed_dec_noise_h6_22 = simmed_stats.mean.BB(si, 3);
        fits(si).simmed_dec_noise_h6_13 = simmed_stats.mean.BB(si, 4);
    
        fits(si).simmed_spatial_bis_h1_22 = simmed_stats.mean.SB(si, 1);
        fits(si).simmed_spatial_bis_h1_13 = simmed_stats.mean.SB(si, 2);
        fits(si).simmed_spatial_bis_h6_22  = simmed_stats.mean.SB(si, 3);
        fits(si).simmed_spatial_bis_h6_13 = simmed_stats.mean.SB(si, 4);
    
        fits(si).simmed_alpha_start = simmed_stats.mean.alpha_start(si);
        fits(si).simmed_alpha_inf = simmed_stats.mean.alpha_inf(si);

     end


    fits = struct2table(fits);
end