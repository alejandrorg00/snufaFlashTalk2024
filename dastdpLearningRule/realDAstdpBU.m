clear; close; clc
rng(100);

M = 100;                % number of synapses per neuron
D = 1;                  % maximal conduction delay
% excitatory neurons    % inhibitory neurons      % total number
Ne = 800;               Ni = 200;                 N = Ne + Ni;
a = [0.02 * ones(Ne,1); 0.1 * ones(Ni,1)];
d = [2 * ones(Ne,1);    2 * ones(Ni,1)];
sm = 4;                 % maximal synaptic strength

post = ceil([N * rand(Ne,M); Ne * rand(Ni,M)]);
s = [ones(Ne,M); -ones(Ni,M)];        % synaptic weights
sd = zeros(N,M);                      % their derivatives
for i = 1:N
    if i <= Ne
        for j = 1:D
            delays{i,j} = M/D * (j - 1) + (1:M/D);
        end
    else
        delays{i,1} = 1:M;
    end
    pre{i} = find(post == i & s > 0);  % pre excitatory neurons
    aux{i} = N * (D - 1 - ceil(ceil(pre{i} / N) / (M/D))) + 1 + mod(pre{i} - 1,N);
end
STDP = zeros(N,1001 + D);
v = -65 * ones(N,1);                   % initial values
u = 0.2 .* v;                          % initial values
firings = [-D 0];                      % spike timings (initialize with dummy spike)

%---------------
% new stuff related to DA-STDP
T = 23;          % the duration of experiment
DA = 0;          % level of dopamine above the baseline
rew = [];

n1 = 1;          % presynaptic neuron
syn = 1;         % the synapse number to the postsynaptic neuron
n2 = post(n1,syn); % postsynaptic neuron
s(n1,syn) = 0;   % start with 0 value

interval = 20;   % the coincidence interval for n1 and n2
n1f = -100;      % the last spike of n1
n2f = [];        % the last spike of n2
shist = zeros(1000 * T, 2);
DAv = zeros(1000*T,1); % DA release

%--------------

% Collect all firings in a separate variable
all_firings = []; % Initialize all_firings to store all spikes

for sec = 1:T                          % simulation of T seconds
    for t = 1:1000                     % simulation of 1 sec
        I = 13 * (rand(N,1) - 0.5);    % random thalamic input
        fired = find(v >= 30);         % indices of fired neurons
        v(fired) = -50;
        u(fired) = u(fired) + d(fired);
        STDP(fired, t + D) = 0.1;
        for kf = 1:length(fired)
            sd(pre{fired(kf)}) = sd(pre{fired(kf)}) + STDP(N * t + aux{fired(kf)});
        end
        firings = [firings; t * ones(length(fired),1), fired];
        % Collect all firings with absolute time
        time = (sec - 1) * 1000 + t; % absolute time in ms
        all_firings = [all_firings; time * ones(length(fired),1), fired];
        k = size(firings,1);
        while firings(k,1) > t - D
            del = delays{firings(k,2), t - firings(k,1) + 1};
            ind = post(firings(k,2), del);
            I(ind) = I(ind) + s(firings(k,2), del)';
            sd(firings(k,2), del) = sd(firings(k,2), del) - 1.5 * STDP(ind, t + D)';
            k = k - 1;
            if k == 0 % To prevent k from becoming zero
                break;
            end
        end
        v = v + 0.5 * ((0.04 * v + 5) .* v + 140 - u + I); % for numerical stability time
        v = v + 0.5 * ((0.04 * v + 5) .* v + 140 - u + I); % step is 0.5 ms
        u = u + a .* (0.2 * v - u);           % recovery variable
        STDP(:, t + D + 1) = 0.95 * STDP(:, t + D); % tau = 20 ms

        DA = DA * 0.995;
        if (mod(t,10) == 0)
            s(1:Ne,:) = max(0, min(sm, s(1:Ne,:) + (0.002 + DA) * sd(1:Ne,:)));
            sd = 0.99 * sd;
        end
        if any(fired == n1)
            n1f = [n1f, sec * 1000 + t];
        end
        if any(fired == n2)
            n2f = [n2f, sec * 1000 + t];
            if (sec * 1000 + t - n1f(end) < interval) && (n2f(end) > n1f(end))
                rew = [rew, sec * 1000 + t + 1000 + ceil(2000 * rand)];
            end
        end
        if any(rew == sec * 1000 + t)
            DA = DA + 0.5;
        end
        DAv((sec - 1) * 1000 + t) = DA;
        shist((sec - 1) * 1000 + t, :) = [s(n1,syn), sd(n1,syn)];
    end
    % Update STDP and firings for the next second
    STDP(:,1:D+1) = STDP(:,1001:1001 + D);
    ind = find(firings(:,1) > 1001 - D);
    firings = [-D 0; firings(ind,1) - 1000, firings(ind,2)];
end
%%
% ---- plot at the end ----
figure;
subplot(4,1,1);
hold on;

% Map neurons n1 and n2 to close y-values to reduce vertical spacing
y_n1 = 1;
y_n2 = 2;

% Get firing times of n1 and n2
times_n1 = all_firings(all_firings(:,2) == n1,1);
times_n2 = all_firings(all_firings(:,2) == n2,1);

% Initialize colors for each spike as gray
colors_n1 = repmat([0.5, 0.5, 0.5], length(times_n1), 1);
colors_n2 = repmat([0.5, 0.5, 0.5], length(times_n2), 1);

% Define the green color for coincidence events
color_green = [51, 204, 51] / 255; % Same green as specified

% Identify coincidence events and change colors to green
for i = 1:length(times_n1)
    t_n1 = times_n1(i);
    % Find indices of n2 spikes that occur after n1 spike but within the interval
    idx_n2 = find(times_n2 > t_n1 & times_n2 <= t_n1 + interval);
    if ~isempty(idx_n2)
        % Change the color of this n1 spike to green
        colors_n1(i,:) = color_green;
        % Change the color of the corresponding n2 spikes to green
        colors_n2(idx_n2,:) = repmat(color_green, length(idx_n2), 1);
    end
end

% Plot vertical bars ('|') for firings of n1
for i = 1:length(times_n1)
    time = times_n1(i);
    color = colors_n1(i,:);
    line([time, time], [y_n1 - 0.4, y_n1 + 0.4], 'Color', color, 'LineWidth', 1);
end

% Plot vertical bars ('|') for firings of n2
for i = 1:length(times_n2)
    time = times_n2(i);
    color = colors_n2(i,:);
    line([time, time], [y_n2 - 0.4, y_n2 + 0.4], 'Color', color, 'LineWidth', 1);
end

% Add labels 'pre' and 'post' next to the neurons
text(-50, y_n1, 'pre', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
text(-50, y_n2, 'post', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

hold off;

% Adjust axis limits to reduce the distance between neurons
axis([0, T * 1000, 0.5, 2.5]); % Y-axis from 0.5 to 2.5 to compress spacing
axis off; % Remove axes

% Add a scale bar indicating 1 s
x0 = 0;
y0 = 2.5;
line([x0, x0 + 1000], [y0, y0], 'Color', 'k', 'LineWidth', 3);
text(x0 + 10, y0 + 0.75, '1 s', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top'); % Label the scale bar

% Plot Synaptic Strength
subplot(4,1,4);
color_synaptic = [5, 20, 53] / 255; % Specified color for synaptic strength
plot(0.001 * (1:(T * 1000)), shist(1:T * 1000,1), 'Color', color_synaptic, 'LineWidth', 1.5);
axis tight; % Adjust axis to fit data
ylim([0 1]);
axis off;   % Remove axes

% Add a scale bar indicating 0.2 mV
hold on;
x0 = T;
y0 = 0.8;            
scale_bar_height = 0.2;                % Height of the scale bar (1 mV)
% Plot the vertical scale bar
plot([x0 x0], [y0 y0 + scale_bar_height], 'k', 'LineWidth', 2);
text(x0 + 0.2, y0 + scale_bar_height / 2, '0.2 mV', 'VerticalAlignment', 'middle', 'FontSize', 10);
hold off;

% Add label 'Synaptic strength' at the beginning of the plot
ylimits = ylim;
x_text = 0.001 * 1; % Beginning of the x-axis
y_text = ylimits(2); % Top of the y-axis
text(x_text, y_text, 'Synaptic strength', 'Color', color_synaptic, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Plot Eligibility Trace
subplot(4,1,3);
color_eligibility = [0.5, 0, 0.5]; % Purple color for eligibility trace
plot(0.001 * (1:(T * 1000)), shist(1:T * 1000,2), 'Color', color_eligibility, 'LineWidth', 1.5);
axis tight;
axis off;

% Add label 'Eligibility trace' at the beginning of the plot
ylimits = ylim;
x_text = 0.001 * 1; % Beginning of the x-axis
y_text = ylimits(2); % Top of the y-axis
text(x_text, y_text, 'Eligibility trace', 'Color', color_eligibility, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Plot Dopamine Level
subplot(4,1,2);
color_dopamine = [51, 204, 51] / 255; % Same green color as used for spikes
plot(0.001 * (1:(T * 1000)), DAv(1:T * 1000), 'Color', color_dopamine, 'LineWidth', 1.5);
axis tight;
axis off;

% Add label 'Dopamine' at the beginning of the plot
ylimits = ylim;
x_text = 0.001 * 1; % Beginning of the x-axis
y_text = ylimits(2); % Top of the y-axis
text(x_text, y_text, 'Dopamine', 'Color', color_dopamine, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

