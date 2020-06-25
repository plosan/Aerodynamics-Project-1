function [NACA, f, p, t, alpha, xh, eta, N] = userInput()

% NACA 4-Series Airfoil
prompt = "Enter the NACA 4-Series Airfoil (e.g. 4412): ";
NACA = -1;
while NACA == -1    
    NACA = str2double(input(prompt, 's'));
    if isnan(NACA)
        fprintf("Non-valid input for NACA Airfoil\n");
        NACA = -1;
    else
        if (NACA < 0000) || (NACA > 9999)
            fprintf("NACA Airfoil out of range\n");
            NACA = -1;
        end
    end
end

f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)
t = mod(NACA, 100)/100;             % Thickness (percent of chord)
if t > 15
    fprintf("NACA %d has a relative thickness of %d%%. The results using TAT may be inaccurate", NACA, t);
end

% Angle of attack
prompt = "Enter the angle of attack of analysis (degrees): ";
alpha = -90;
while alpha == -90
    alpha = str2double(input(prompt, 's'));
    if isnan(alpha)
        fprintf("Non-valid input for angle of attack\n");
        NACA = -1;
    else
        if (alpha < -15) || (alpha > 15)
            fprintf("Angle of attack out of linear range\n");
            alpha = -90;
        end
    end
end

% Flap hinge location
prompt = "Enter the flap hinge location (0 <= xh <= 1): ";
xh = -1;
while xh == -1
    xh = str2double(input(prompt, 's'));
    if isnan(xh)
        fprintf("Non-valid input for flap hinge location\n");
        xh = -1;
    else
        if (xh < 0) || (xh > 1)
            fprintf("Flap hinge location out of range\n");
            xh = -1;
        end
    end
end

% Flap deflection angle
prompt = "Enter the flap deflection angle (degrees): ";
eta = -90;
while eta == -90
    eta = str2double(input(prompt, 's'));
    if isnan(xh)
        fprintf("Non-valid input for flap deflection angle\n");
        eta = -1;
    else
        if (eta < -15) || (eta > 15) %% REVISAR
            fprintf("Flap deflection angle out of range\n");
            eta = -90;
        end
    end
end

% Number of stages
prompt = "Enter the number of stages (0 < N < 201): ";
N = -1;
while N == -1
    N = str2double(input(prompt, 's'));
    if isnan(xh)
        fprintf("Non-valid input for number of stages\n");
        N = -1;
    else
        N = round(N);
        if (N < 1) || (N > 200)
            fprintf("Number of stages out of range\n");
            N = -1;
        end
    end
end

end