function [sc_param, dr_sc_param, P, dPdr] = MARGO_param(r)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP constant parameters

% MARGO mission adimensional propulsive parameters, and gradients, given adimensional distance

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    MU=22.6;                                    % m0 [kg]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    g0=9.80665e-3;      % [km/s^2]

    a=[-0.7253 0.02481 0 0 0];
    b=[2652 -18.123 0.3887 -0.00174 0];

    Tmax=@(P) sum(P.^(0:4).*a,2); % [mN]
    Isp=@(P) sum(P.^(0:4).*b,2);  % [s]
    Pin=@(r) 120;

    sc_param=zeros(2,length(r));
    dr_sc_param=zeros(2,length(r));

    sc_param(1,:)=TU^2/(MU*LU)*Tmax(Pin(r))*1e-6;  % T [-]
    sc_param(2,:)=TU/LU*Isp(Pin(r))*g0;            % c [-]

    P=120*ones(size(r));

    dPdr=0;
    dr_sc_param(1,:)=0;
    dr_sc_param(2,:)=0;

end