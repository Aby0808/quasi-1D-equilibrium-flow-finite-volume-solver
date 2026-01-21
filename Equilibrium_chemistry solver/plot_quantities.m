%% ============================
%  Plot HEGEL quasi-1D nozzle .dat
%  ============================

clear; clc; close all;

% --------- USER INPUT ----------
fname = 'F:\study\UIUC_pg\Thesis_MS\04_output\quasi_1d_nozzle\test\flowfield_EU_NLTE.dat';   % <-- change to your file name
% -------------------------------

%% 1) Read file as text, remove headers, fix wrapped lines ending with '\'
raw = fileread(fname);

% Split into lines
lines = regexp(raw, '\r\n|\n|\r', 'split');

% Remove comment/header lines starting with '#'
lines = lines(~startsWith(strtrim(lines), '#'));
lines = lines(~cellfun(@isempty, strtrim(lines)));

% HEGEL sometimes wraps long lines with a trailing '\'
% Merge lines that end with '\' with the following line
merged = {};
i = 1;
while i <= numel(lines)
    L = strtrim(lines{i});
    while ~isempty(L) && L(end) == '\'
        L = strtrim(L(1:end-1)); % remove trailing '\'
        if i+1 <= numel(lines)
            L = [L, ' ', strtrim(lines{i+1})]; %#ok<AGROW>
            i = i + 1;
        else
            break;
        end
    end
    merged{end+1,1} = L; %#ok<SAGROW>
    i = i + 1;
end

% Join back to one big numeric string and parse
numStr = strjoin(merged, newline);
vals = sscanf(numStr, '%f');

%% 2) Define columns (based on your header)
% Header says:
% x  X_em  X_N  X_O  X_N2  X_NO  X_O2  X_N2p  X_NOp  X_Np  X_O2p  X_Op  Th  Tve  rho  p  H  Mf  u
ncol = 19;

if mod(numel(vals), ncol) ~= 0
    error(['Parsed %d numbers, which is not divisible by %d columns.\n' ...
           'This usually means the file has extra tokens/formatting issues.'], numel(vals), ncol);
end

data = reshape(vals, ncol, []).';   % rows = x stations

% Assign columns
x   = data(:,1);
Xem = data(:,2);
XN  = data(:,3);
XO  = data(:,4);
XN2 = data(:,5);
XNO = data(:,6);
XO2 = data(:,7);
XN2p= data(:,8);
XNOp= data(:,9);
XNp = data(:,10);
XO2p= data(:,11);
XOp = data(:,12);
Th  = data(:,13);
Tve = data(:,14);
rho = data(:,15);
p   = data(:,16);
H   = data(:,17);
Mf  = data(:,18);
u   = data(:,19);

% Sort by x just in case
[x, is] = sort(x);
Xem=Xem(is); XN=XN(is); XO=XO(is); XN2=XN2(is); XNO=XNO(is); XO2=XO2(is);
XN2p=XN2p(is); XNOp=XNOp(is); XNp=XNp(is); XO2p=XO2p(is); XOp=XOp(is);
Th=Th(is); Tve=Tve(is); rho=rho(is); p=p(is); H=H(is); Mf=Mf(is); u=u(is);

%% 3) Plots
set(groot, 'defaultLineLineWidth', 1.8);
set(groot, 'defaultAxesFontSize', 12);

% ---- Temperatures
figure('Name','Temperatures');
plot(x, Th, '-'); hold on;
plot(x, Tve, '--');
grid on; xlabel('x (m)'); ylabel('T (K)');
legend('T_h','T_{ve}','Location','best');
title('Temperatures vs x');

% ---- Flow variables
figure('Name','Flow variables');
yyaxis left
plot(x, p, '-'); ylabel('p (Pa)');
yyaxis right
plot(x, u, '--'); ylabel('u (m/s)');
grid on; xlabel('x (m)');
title('Pressure and velocity vs x');
legend('p','u','Location','best');

%% ---- Density
figure('Name','Density');
plot(x, rho, '-', 'LineWidth', 1.8);
grid on;
xlabel('x (m)');
ylabel('\rho (kg/m^3)');
title('Density vs x');

%% ---- Enthalpy
figure('Name','Enthalpy');
plot(x, H, '--', 'LineWidth', 1.8);
grid on;
xlabel('x (m)');
ylabel('H (J/kg)');
title('Enthalpy vs x');

%% ---- Mach number
figure('Name','Mach number');
plot(x, Mf, '-.', 'LineWidth', 1.8);
grid on;
xlabel('x (m)');
ylabel('M_f');
title('Mach number vs x');


% ---- Species (pick the ones you care about)
figure('Name','Species mass fractions');
plot(x, XN2, '-'); hold on;
plot(x, XN, '--');
plot(x, XO, '--');
plot(x, XNO, '-.');
plot(x, XO2, '-.');
grid on; xlabel('x (m)'); ylabel('X');
legend('X_{N2}','X_N','X_O','X_{NO}','X_{O2}','Location','best');
title('Neutral species vs x');

figure('Name','Ions / electrons');
plot(x, Xem, '-'); hold on;
plot(x, XN2p, '--');
plot(x, XNOp, '--');
plot(x, XNp, '-.');
plot(x, XO2p, '-.');
plot(x, XOp, '-.');
grid on; xlabel('x (m)'); ylabel('X');
legend('X_{em}','X_{N2+}','X_{NO+}','X_{N+}','X_{O2+}','X_{O+}','Location','best');
title('Charged species vs x');

