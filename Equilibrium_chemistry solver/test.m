%% ============================================
%  Quick sanity plots for equilibrium tables
% ============================================

T = table.T;      % NT
P = table.P;      % NP

% Pick representative indices
iP = round(length(P)/2);   % mid pressure
iT = round(length(T)/2);   % mid temperature

%% ---- gamma vs T at fixed P ----
figure;
plot(T, table.gamma(:, iP), 'LineWidth', 1.5);
grid on;
xlabel('Temperature (K)');
ylabel('\gamma');
title(sprintf('\\gamma vs T at P = %.2e Pa', P(iP)));

%% ---- Cp vs T at fixed P ----
figure;
plot(T, table.Cp(:, iP), 'LineWidth', 1.5);
grid on;
xlabel('Temperature (K)');
ylabel('C_p (J/kg-K)');
title(sprintf('C_p vs T at P = %.2e Pa', P(iP)));

%% ---- gamma vs P at fixed T ----
figure;
semilogx(P, table.gamma(iT, :), 'LineWidth', 1.5);
grid on;
xlabel('Pressure (Pa)');
ylabel('\gamma');
title(sprintf('\\gamma vs P at T = %.0f K', T(iT)));

%% ---- Cp vs P at fixed T ----
figure;
semilogx(P, table.Cp(iT, :), 'LineWidth', 1.5);
grid on;
xlabel('Pressure (Pa)');
ylabel('C_p (J/kg-K)');
title(sprintf('C_p vs P at T = %.0f K', T(iT)));


%% ============================================
%  2D contour plots in (T,P)
% ============================================

[TG, PG] = meshgrid(T, P);   % note transpose needed for plotting

%% ---- gamma(T,P) ----
figure;
contourf(T, P, table.gamma.', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title('\gamma(T,P)');
grid on;

%% ---- Cp(T,P) ----
figure;
contourf(T, P, table.Cp.', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title('C_p(T,P) [J/kg-K]');
grid on;

%% ---- Mmix(T,P) ----
figure;
contourf(T, P, table.M.', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title('Mean molecular weight (kg/mol)');
grid on;

%% ============================================
%  Species mole fraction check
% ============================================

isp = 1;  % choose species index (e.g. N2)
figure;
contourf(T, P, squeeze(table.X(isp,:,:)).', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title(sprintf('N2 mole fraction'));
grid on;

isp = 2;  % choose species index (e.g. N2)
figure;
contourf(T, P, squeeze(table.X(isp,:,:)).', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title(sprintf('N mole fraction'));
grid on;

isp = 3;  % choose species index (e.g. N2)
figure;
contourf(T, P, squeeze(table.X(isp,:,:)).', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title(sprintf('O2 mole fraction'));
grid on;

isp = 4;  % choose species index (e.g. N2)
figure;
contourf(T, P, squeeze(table.X(isp,:,:)).', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title(sprintf('O mole fraction'));
grid on;

isp = 5;  % choose species index (e.g. N2)
figure;
contourf(T, P, squeeze(table.X(isp,:,:)).', 30, 'LineColor', 'none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Temperature (K)');
ylabel('Pressure (Pa)');
title(sprintf('NO mole fraction'));
grid on;