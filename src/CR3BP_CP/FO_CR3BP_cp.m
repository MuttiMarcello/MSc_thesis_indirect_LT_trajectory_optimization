function [dydPhi] = FO_CR3BP_cp(~, z, prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% implicit time optimal energy-to-fuel optimal problem dynamics: Euler-Lagrange
% equations + STM dynamics

    % propulsive parameters
    U=prob.param;

    y=z(1:14); % complete state
    Phi=reshape(z(15:210),[14,14]); % STM

    [dy,A]=FO_SEP_dz(y,U);
        
    dydPhi=zeros(210,1);

    dydPhi(1:14)=dy;
    dydPhi(15:210)=reshape(A*Phi,[14*14,1]);

end

function [dy,A]=FO_SEP_dz(y,U)

        [dy,A]=u_on_dz(y,U);    % u=1

end