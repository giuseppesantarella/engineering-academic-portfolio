function results = final_controller_1deputy2(currentx, currentr, Ts, t, A_vec, B_vec)
% MPC LPV con YALMIP ottimizzando solo su u_y e u_z
% u_x = 0 impostato a monte

persistent ctrl nx nu ny N theta_max u_min u_max N_c a_max A_max pastu Deltau_min Deltau_max ctrl2 ctrl3

%% Parametri statici
if isempty(nx)
    nx = 7;       % dimensione stato
    nu = 2;       % SOLO u_y, u_z
    ny = nx;
    N  = 5;
    N_c = 5;
    theta_max = 45*pi/180;
    pastu = zeros(3,1);
end

% Trasformazione input
B_vec = B_vec.*1e-6; % Trasformato in secondi/millimetro

%% Inizializzazione controller
if isempty(ctrl)
    max_thrust = 0.65; % [milliNewton]
    mass_d = 20;        % [kg]
    u_max = max_thrust / mass_d; % [mm/s^2]
    min_thrust = 0.35;
    u_min = min_thrust / mass_d; % [mm/s^2]
    Deltau_max=0.51*u_min;
    Deltau_min=-Deltau_max;

    A_max = 1.95*u_min;
    a_max = 0.95*A_max;

    % Variabili decisionali
    X = sdpvar(nx, N+1);
    X_u=sdpvar(nu,N+1);
    Ured = sdpvar(nu, N_c);   % [u_y; u_z] decisione
    DeltaUred=sdpvar(nu,N_c);
    x0 = sdpvar(nx,1);
    r  = sdpvar(nx,1);
    u_past=sdpvar(3,1);

    % Parametri numerici (input) per YALMIP
    Adp = sdpvar(nx,nx);
    Bdp_eff = sdpvar(nx,nu);
    Adp_g=[Adp,Bdp_eff;
        zeros(nu,nx),eye(nu)];
    Bdp_eff_g=[Bdp_eff;
        eye(nu)]; % matrice globale degli input

    % Funzione costo di default
  
    Q = eye(nx);
    e_x = [1e-5; 1e-6; 1e-6; 1e-6; 1e-6; 1e-6]; % qui le unità degli elementi di e_x sono kilometri
    for j = 1:6
        Q(j,j) = Q(j,j)./(e_x(j).^2);
    end
    % Modifica euristica dei parametri
    Q(2,2)=Q(2,2).*1700;
    Q(3,3)=Q(3,3).*20;
    Q(4,4)=Q(4,4).*30;
    Q(5,5)=Q(5,5).*10;

    R = eye(nu);
    R(1,1)=R(1,1)./(u_max^2); % le unità di misura di u_max qui sono mm/s^2 
    R(2,2)=R(2,2)./(u_max^2);
    Q_g=[Q, zeros(nx,nu);
        zeros(nu,nx),R];

    R_g = eye(nu);
    R_g(1,1)=R_g(1,1)./(Deltau_max^2); % le unità di misura di Deltau_max qui sono mm/s^2 
    R_g(2,2)=R_g(2,2)./(Deltau_max^2);

    constraints = [];
    objective = 0;

    constraints=[constraints,X_u(:,1)==u_past(2:3)];
    for k = 1:N
        constraints=[constraints, Deltau_min<=DeltaUred(:,k), DeltaUred(:,k)<=Deltau_max];
        if k>=2
            % constraints = [constraints, 0 <= X_u(1,k), X_u(1,k) <= u_max];
            constraints = [constraints, -u_max <= X_u(1,k), X_u(1,k) <= u_max];
            constraints = [constraints, -u_max <= X_u(2,k), X_u(2,k) <= u_max];
        end

        objective = objective + ([r;zeros(nu,1)]-[X(:,k);X_u(:,k)])'*Q_g*([r;zeros(nu,1)]-[X(:,k);X_u(:,k)]) + DeltaUred(:,k)'*R_g*DeltaUred(:,k);
        constraints = [constraints, [X(:,k+1);X_u(:,k+1)] == Adp_g*[X(:,k);X_u(:,k)] + Bdp_eff_g*DeltaUred(:,k)];
    end
    constraints = [constraints, -u_max <= X_u(:,N+1), X_u(:,N+1) <= u_max];
    objective = objective + ([r;zeros(nu,1)]-[X(:,N+1);X_u(:,N+1)])'*Q_g*([r;zeros(nu,1)]-[X(:,N+1);X_u(:,N+1)]);
    constraints = [constraints, X(:,1) == x0];

    ops = sdpsettings('solver','mosek','verbose',1,'debug',1); % solutore
    % in uso attualmente

    % ops = sdpsettings('solver','osqp','verbose',2,'debug',1);

    % ops = sdpsettings('solver','gurobi','verbose',2,'debug',1);

    % ops = sdpsettings('solver','mosek','verbose',2,'debug',1, ...
    %     'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-12, ...
    %     'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS',1e-12);

    %
    % ops = sdpsettings('solver','osqp','verbose',2,'debug',1, ...
    %     'osqp.eps_abs', 1e-8, 'osqp.eps_rel', 1e-7, 'osqp.max_iter', 10000);

    % Creo optimizer
    % ctrl = optimizer(constraints, objective, ops, {x0, r, Adp, Bdp_eff}, Ured);
    ctrl = optimizer(constraints, objective, ops, {x0, r, Adp, Bdp_eff, u_past}, X_u);
    ctrl2=optimizer(constraints, objective, ops, {x0, r, Adp, Bdp_eff, u_past}, X_u);
    ctrl3=optimizer(constraints, objective, ops, {x0, r, Adp, Bdp_eff, u_past}, DeltaUred);
end

%% Ricostruisco Ad e Bd dal vettore
A = column_vector2matrix(A_vec, nx);
B = column_vector2matrix(B_vec, 3);
B_eff = B(:,2:3); % solo u_y, u_z
C = eye(nx);
D = zeros(ny,3);

% Discretizzazione
sysc_ss = ss(A,B,C,D);
sysd_ss = c2d(sysc_ss, Ts);
Ad = sysd_ss.A;
Bd = sysd_ss.B;
Bd_eff = sysd_ss.B(:,2:3); % è identico rispetto a fare Bd = sysd_ss.B; Bd_eff = Bd(:,2:3);
% format shortE
% disp(Bd_eff)
if Bd(:,2:3) ~= Bd_eff
    disp('le matrici B sono diverse')
end
Cd = C;

%% Calcolo uscita stimata
y_est = Cd*currentx;

%% Calcolo MPC
% [ured_seq, errorcode] = ctrl({currentx, currentr, Ad, Bd_eff, pastu});
[X_u_seq, errorcode] = ctrl({currentx, currentr, Ad, Bd_eff, pastu});
% [X_u, errorcode2] = ctrl2({currentx, currentr, Ad, Bd_eff, pastu});
[Deltau_seq, errorcode3] = ctrl3({currentx, currentr, Ad, Bd_eff, pastu});
if errorcode ~= 0
    yalmiperror(errorcode)
    error('MPC:SolverFailure', 'YALMIP failed at t = %g, problem = %d.', t, errorcode);
end
% uout = [0; ured_seq(:,1)];
uout = [0; X_u_seq(:,2)];

% pastu
% ured_seq
% X_u_seq

% seq2=[];
% pastu
% sum2=pastu(2:3)+Deltaured_seq(:,1);
% seq2(:,1)=sum2;
% for k=2:N
%     sum2=sum2+Deltaured_seq(:,k);
%     seq2(:,k)=sum2;
% end
% seq2

% seq3=[];
% for k=1:N
%     seq3(:,k)=X_u(:,k+1);
% end
% seq3

% Deltau_seq
tol = 1e-6; % ad esempio tolleranza numerica
%
% if max(abs(ured_seq(:,1) - (pastu(2:3) + Deltau_seq(:,1)))) > tol
%     disp('qualcosa non va nel controllore al passo 1')
% end
%
% for k = 2:N
%     if max(abs(ured_seq(:,k) - (ured_seq(:,k-1) + Deltau_seq(:,k)))) > tol
%         disp('qualcosa non va nel controllore al passo %d',k)
%     end
% end

if max(abs(X_u_seq(:,2) - (pastu(2:3) + Deltau_seq(:,1)))) > tol
    disp('qualcosa non va nel controllore al passo 1')
end

for k = 2:N
    if max(abs(X_u_seq(:,k+1) - (X_u_seq(:,k) + Deltau_seq(:,k)))) > tol
        disp('qualcosa non va nel controllore al passo %d',k)
    end
end

%% CHECK VINCOLI POST CALCOLO DI OTTIMIZZAZIONE
tol = 2e-3; % [mm/s^2] questa è secondo me la quantità minima oltre la
% violazione dei vincoli diventa eccessiva, considerando che u_max è pari a
% circa 3.2e-2 mm/s^2.

violazioni = strings(0,1);   % string array vuoto

for k = 1:N
    % u = [0; ured_seq(:,k)];  % ricostruisco u completo (3x1)
    u = [0; X_u_seq(:,k+1)];

    % 1) moduli delle componenti di u

    viol=abs(u(2)) - u_max;
    if viol > tol
        violazioni(end+1) = sprintf('Violazione modulo 2^ componente di u a step %d della quantità %d mm/s^2', k, viol);
    end

    viol=abs(u(3)) - u_max;
    if viol > tol
        violazioni(end+1) = sprintf('Violazione modulo 3^ componente di u a step %d della quantità %d mm/s^2', k, viol);
    end

    % 2) u_y >= 0
    % viol_u_y = u(2) + tol;
    % if viol_u_y < 0
    %     violazioni(end+1) = sprintf('Violazione u_y>=0 a step %d: tolleranza sforata di %d mm/s^2', k, viol_u_y);
    % end

    % % 3) angolo (applicato solo in caso di station keeping)
    % viol = abs(u(3)) - abs(u(2));
    % if abs(u(3)) > abs(u(2))*tan(theta_max) + tol
    %     violazioni(end+1) = sprintf('Violazione angolo a step %d pari a %d', k, viol);
    % end

    % if k == 1
    %     % fprintf('pastu è %d, %d, %d [mm/s^2]. \n', pastu(1), pastu(2), pastu(3));
    %
    %     % viol = abs(ured_seq(1,k)-pastu(2)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew u_y tra 1° step e pastu di %d mm/s^2', viol);
    %     % end
    %     %
    %     % viol = abs(ured_seq(2,k)-pastu(3)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew u_z tra 1° step e pastu di %d mm/s^2', viol);
    %     % end
    %
    %     % viol = abs(ured_seq(1,k)) + abs(pastu(3)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew incrociati tra 1° step e pastu di %d mm/s^2', viol);
    %     % end
    %     %
    %     % viol = abs(ured_seq(2,k)) + abs(pastu(2)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew incrociati tra 1° step e pastu di %d mm/s^2', viol);
    %     % end
    %
    % end

    % 4) slew rate (solo se non è l’ultimo)
    % if k <= N_c-1
    %
    %     % viol = abs(ured_seq(1,k+1)-ured_seq(1,k)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew u_y tra step %d e %d della quantità %d mm/s^2', k,k+1, viol);
    %     % end
    %     %
    %     % viol = abs(ured_seq(2,k+1)-ured_seq(2,k)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew u_z tra step %d e %d della quantità %d mm/s^2', k,k+1, viol);
    %     % end
    %
    %     % viol = abs(ured_seq(1,k+1)-ured_seq(2,k)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew incrociati tra step %d e %d della quantità %d mm/s^2', k,k+1, viol);
    %     % end
    %     %
    %     % viol = abs(ured_seq(2,k+1)-ured_seq(1,k)) - A_max;
    %     % if viol > tol
    %     %     violazioni(end+1) = sprintf('Violazione slew incrociati tra step %d e %d della quantità %d mm/s^2', k,k+1,viol);
    %     % end
    %
    % end
    % if k == N_c
    %     viol = abs(ured_seq(1,k)-ured_seq(2,k)) - A_max;
    %     if viol > tol
    %         violazioni(end+1) = sprintf('Violazione slew incrociati allo step %d della quantità %d mm/s^2', k, viol);
    %     end
    % end
end

if ~isempty(violazioni)
    warning('Ad istante %d s sono state trovate le seguenti violazioni dei vincoli:', t);
    for i = 1:numel(violazioni)
        fprintf('%s\n', violazioni(i));
    end
    % else
    %     disp('Tutti i vincoli rispettati');
end


%% Saturazione per il lower bound
% if 0.5*u_min<abs(uout(2)) && abs(uout(2))<=u_min
%     uout(2)=u_min;
% end
% if abs(uout(2))<=0.5*u_min
%     uout(2)=0;
% end
%
% if 0.5*u_min<abs(uout(3)) && abs(uout(3))<=u_min
%     uout(3)=u_min;
% end
% if abs(uout(3))<=0.5*u_min
%     uout(3)=0;
% end

%% Output finale
pastu = uout;
uout = uout*1e-6; % mm/s^2 -> km/s^2
results = [uout; y_est];

end


