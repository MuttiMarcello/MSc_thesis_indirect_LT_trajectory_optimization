function [tout,xout] = TO_ode87(prob,tspan,z0)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% Simplified integration scheme for TO problems
% used build in ode as no switches are anticipated

    odopt = odeset(RelTol=5e-14,AbsTol=1e-14);
    [tout,xout]=ode113(@(t,y) TO_CR3BP_cp(t,y,prob),tspan,z0,odopt);

end