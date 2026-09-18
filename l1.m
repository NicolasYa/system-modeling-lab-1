clear; clc; close all;

%% параметры
k  = 1;   l  = 10;  m  = 1;   n  = 10;
kt = 100; b  = 22000; i1 = 10; i2 = 1;
s  = 200; V  = 800;  T  = 14;
x0 = [0.2; 0.2; 0; 0; 0];          % [x1; x2; x3; delta; x5]
delta_max = 0.5;

h = input('Введите шаг интегрирования h: ');
if isempty(h) || h <= 0
    h = 0.01;
    fprintf('Некорректный ввод. Установлен шаг по умолчанию: h = %g\n', h);
end

fprintf('Решение с заданным шагом h = %.6f\n', h);

[t, X, x5_T] = euler_model2(h, T, x0, k,l,m,n,kt,b,i1,i2,s,V,delta_max);

% эталон с половинным шагом
[~, ~, x5_T2] = euler_model2(h/2, T, x0, k,l,m,n,kt,b,i1,i2,s,V,delta_max);
delta_rel = abs(x5_T2 - x5_T)/abs(x5_T2)*100;

fprintf('x5(T) = %.10f\n', x5_T);
fprintf('Относительная погрешность δ = %.4f %%\n\n', delta_rel);

% график 1
figure('Name','Решение системы с заданным h', 'NumberTitle','off');
subplot(5,1,1);
plot(t, X(1,:), 'b', 'LineWidth', 1.3);
ylabel('x_1'); grid on;
subplot(5,1,2);
plot(t, X(2,:), 'r', 'LineWidth', 1.3);
ylabel('x_2'); grid on;
subplot(5,1,3);
plot(t, X(3,:), 'g', 'LineWidth', 1.3);
ylabel('x_3'); grid on;
subplot(5,1,4);
plot(t, X(4,:), 'm', 'LineWidth', 1.3);
ylabel('\delta'); grid on;
subplot(5,1,5);
plot(t, X(5,:), 'k', 'LineWidth', 1.3);
ylabel('x_5'); xlabel('t, с'); grid on;
sgtitle(sprintf('h = %.4f | x_5(T) = %.6f | δ = %.3f%%', ...
              h, x5_T, delta_rel));

%% анализ зависимости точности и трудоёмкости 
fprintf('Анализ зависимости δ и трудоёмкости от шага h \n');

h_list = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001];
err_list   = zeros(size(h_list));
steps_list = zeros(size(h_list));

fprintf('%-12s %-18s %-12s\n', 'h', 'x5(T)', 'δ, %');
fprintf('--------------------------------------------\n');

for i = 1:length(h_list)
    h = h_list(i);
    [~, ~, y]  = euler_model2(h,   T, x0, k,l,m,n,kt,b,i1,i2,s,V,delta_max);
    [~, ~, y2] = euler_model2(h/2, T, x0, k,l,m,n,kt,b,i1,i2,s,V,delta_max);
    
    err_list(i)   = abs(y2 - y)/abs(y2)*100;
    steps_list(i) = ceil(T/h);
    
    fprintf('%-12.6f %-18.8f %-12.4f\n', h, y, err_list(i));
end
fprintf('\n');

% график 2
figure('Name','Зависимость погрешности и трудоёмкости', 'NumberTitle','off');
subplot(2,1,1);
loglog(h_list, err_list, 'b-o', 'LineWidth', 1.6, 'MarkerSize', 7);
grid on; xlabel('h'); ylabel('\delta, %');
title('Относительная погрешность по x_5(T)');
subplot(2,1,2);
loglog(h_list, steps_list, 'r-s', 'LineWidth', 1.6, 'MarkerSize', 7);
grid on; xlabel('h'); ylabel('Число шагов J = T/h');
title('Трудоёмкость (число шагов)');

%% автоматический выбор шага (δ ≤ 1 %)
fprintf('Автоматический выбор шага (требуется δ ≤ 1 %%)\n');
fprintf('%-12s %-18s %-12s\n', 'h', 'x5(T)', 'δ, %');
fprintf('--------------------------------------------\n');

h = 0.1; % начальное значение
max_iter = 25;

for iter = 1:max_iter
    [t, X, y]  = euler_model2(h,   T, x0, k,l,m,n,kt,b,i1,i2,s,V,delta_max);
    [~, ~, y2] = euler_model2(h/2, T, x0, k,l,m,n,kt,b,i1,i2,s,V,delta_max);
    
    delta_rel = abs(y2 - y)/abs(y2)*100;
    
    fprintf('%-12.8f %-18.8f %-12.4f\n', h, y, delta_rel);
    
    if delta_rel <= 1.0
        break;
    end
    h = h / 2;
end

fprintf('\n------------------------------------------------------------\n');
fprintf('Итоговый результат:\n');
fprintf('  Выбранный шаг h          = %.8f\n', h);
fprintf('  x5(T)                    = %.10f\n', y);
fprintf('  Относительная погрешность δ = %.4f %%\n', delta_rel);
fprintf('  Число шагов              = %d\n', ceil(T/h));
fprintf('------------------------------------------------------------\n\n');

% график 3
figure('Name','Итоговое решение с автоматически подобранным h', 'NumberTitle','off');
subplot(5,1,1);
plot(t, X(1,:), 'b', 'LineWidth', 1.3);
ylabel('x_1'); grid on;
subplot(5,1,2);
plot(t, X(2,:), 'r', 'LineWidth', 1.3);
ylabel('x_2'); grid on;
subplot(5,1,3);
plot(t, X(3,:), 'g', 'LineWidth', 1.3);
ylabel('x_3'); grid on;
subplot(5,1,4);
plot(t, X(4,:), 'm', 'LineWidth', 1.3);
ylabel('\delta'); grid on;
subplot(5,1,5);
plot(t, X(5,:), 'k', 'LineWidth', 1.3);
ylabel('x_5'); xlabel('t, с'); grid on;
sgtitle(sprintf('автоматически выбранный h = %.6f | δ = %.3f%%', ...
              h, delta_rel));

%% функция метода Эйлера
function [t_out, X_out, x5_T] = euler_model2(h, T, x0, k,l,m,n,kt,b,i1,i2,s,V,delta_max)
    N = ceil(T/h) + 1;
    t = zeros(1, N);
    X = zeros(5, N);
    X(:,1) = x0;
    t(1) = 0;
    
    % сохраняются точки каждые 0.01 с
    save_step = max(1, round(0.01 / h));
    idx = 1;
    t_out = t(1);
    X_out = X(:,1);
    
    for j = 1:N-1
        tt = t(j);
        xx = X(:,j);
        
        % насыщение δ->x4
        if abs(xx(4)) <= delta_max
            x4 = xx(4);
        else
            x4 = delta_max * sign(xx(4));
        end
        
        % θ
        den = b - V * tt;
        if abs(den) < 1e-12
            theta = 0;
        else
            theta = (10000 - xx(5)) / den;
        end
        
        % правые части
        dx1    = k * (xx(2) - xx(1));
        dx2    = xx(3);
        dx3    = l*(xx(1) - xx(2)) - m*xx(3) + n*x4;
        ddelta = -kt*x4 - i1*xx(2) - i2*xx(3) + s*(theta - xx(2));
        dx5    = V * sin(xx(1));
        
        % шаг Эйлера
        X(:,j+1) = xx + h * [dx1; dx2; dx3; ddelta; dx5];
        t(j+1)   = tt + h;
        
        % сохранение для графика
        if mod(j, save_step) == 0 || j == N-1
            idx = idx + 1;
            t_out(idx) = t(j+1);
            X_out(:,idx) = X(:,j+1);
        end
    end
    x5_T = X(5, end);
end