%% This uses finite differences to solve diffusion
function vals = diffDiffusion(endtime, timestep, interval, spacestep, ...
    initData)
    % Approximates diffusion by finite differences
    % Make large array for values at each point in time and space
    % timestep / (spacestep)^2 is supposed to be less than 1/2
    timesteps = endtime / timestep;
    spacesteps = (interval(2) - interval(1)) / spacestep;
    s = timestep / (spacestep)^2; % parameter for solver
    vals = zeros(timesteps, spacesteps);

    % Value at boundary is zero at each time.
    vals(:,1) = zeros(timesteps,1);
    vals(:, end) = zeros(timesteps,1);

   % Add the initial values.
    vals(1,:) = arrayfun(initData, linspace(interval(1),interval(2), spacesteps));

    for i = 2:timesteps
        for j = 2:(spacesteps -1)
            vals(i, j) = s * (vals(i-1, j+1) + vals(i-1, j-1)) ...
            + (1-2*s) * vals(i-1, j);
        end
    end
end

endtime = 0.5;
interval = [0,pi];
timestep = 0.005;
spacestep = pi/20; % should guarantee stability-ish
initData = @(x) x * (pi - x);

function y = initData2(x)
    if (0 <= x) && (x < 1)
        y = x;
    else
        y = 0;
    end
end

dumbData = @(x) initData2(x);

g = diffDiffusion(endtime, timestep, interval, spacestep, initData);
dg = diffDiffusion(endtime, timestep, interval, spacestep, dumbData);

figure(1)
mesh(g)
figure(2)
mesh(dg)

%% This uses the Fourier method instead

function coeff = fouCoeff(func, freq, dom)
    % dom is the domain, an interval
    % func is the function to analyze, as a handle
    % freq is the frequency at which we take the coefficient
    % returns the boring sine series coefficients
    l = dom(2) - dom(1);
    intfunc = @(x) func(x) .* sin(freq * pi * x / l);
    coeff = 2 / l .* integral(intfunc, dom(1), dom(2));
end

function y = diffKer(coeffs, len, l, x,t)
    y = 0;
    for i = 1:1:len
        y = y + coeffs(i) .* exp(-1*(i * pi / l)^2 * t) ...
            * sin(pi * i .* x / l);
    end
end

function vals = fouDiffusion(endtime, timestep, interval, spacestep, initData, deg)
    length = interval(2)-interval(1);
    timesteps = endtime/ timestep;
    spacesteps = length / spacestep;
    fouCoeffs = zeros(deg, 1);
    for i = 1:1:deg
        fouCoeffs(i, 1) = fouCoeff(initData, i, interval);
    end
    vals = zeros(timesteps, spacesteps);

    % Function that computes the solution at (x,t)
    sol = @(x,t) diffKer(fouCoeffs, deg, length,x,t);
    
    fouVals = zeros(spacesteps, timesteps, 2); % array of time/space vals
    for i = 1:1:spacesteps
        for j = 1:1:timesteps
            fouVals(i , j, 1) = i * length / spacesteps;
            fouVals(i, j, 2) = j * endtime / timesteps;
        end
    end

    for i = 1:1:spacesteps
        for j = 1:1:timesteps
            vals(j,i) = sol(fouVals(i,j,1),fouVals(i,j,2));
        end
    end

    % Value at boundary is zero at each time.
    vals(:,1) = zeros(timesteps,1);
    vals(:, end) = zeros(timesteps,1);

   % Add the initial values.
    vals(1,:) = arrayfun(initData, linspace(interval(1),interval(2), ...
        spacesteps));
end

% same initData as before
endtime = 0.5;
timestep = 0.01;
interval = [0,pi];
spacestep = pi/10;
deg = 7;
initDataNew = @(x) x .* (pi -x);

vals = fouDiffusion(endtime, timestep, interval, spacestep, ...
initDataNew, deg);

figure(3)
mesh(vals)
%%

quotient = timestep/spacestep^2 % we hold the quotient constant


