function prob = struct_assembly(twin,m0,param,x0,xf,isFO)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% problem structure field assembly

    % targ    [str]
    % m0      [1x1, kg]
    % Plim    [1x2, W, (Pmin Pmax)]
    % twin    [1x2, s, (t_wo t_wc)]
    % isFO    [1x1]
    
    %-problem-setup--------------------------------------------------------
    prob.tw=twin;           % [1x2, s]
    prob.m0=m0;             % [1x1, kg]
    prob.param=param;       % [1x3, -] [T, c, mu]
    prob.x0=x0;             % [3x1, -]
    prob.xf=xf;             % [3x1, -]
    prob.isFO=isFO;         % [1x1, -]
    
    if isFO==0
        prob.epsilon=0;
    else
        prob.epsilon=1;    % EFO continuation parameter
        prob.isTO=0;
    end
    
    %-added-by-TO_t0CONT---------------------------------------------------
    prob.t0=[];             % [1x1, s]
    
    %-added-by-ZFP---------------------------------------------------------
    prob.gamma=[];          % [8x1, -]
    prob.jac=[];            % [8x8, -]
    prob.tf=[];             % [1x1, s]
    prob.tf_ad=[];          % [1x1, -]
    prob.y0=[];             % [14x1, -]

    %-added-by-DispRes-----------------------------------------------------
    prob.zz=[];             % [Nx210, -]
    prob.tt=[];             % [Nx1, d]
    prob.tt_ad=[];          % [Nx1, -]
    prob.mf=[];             % [1x1, kg]
    prob.mp=[];             % [1x1, kg]
    prob.S=[];              % [Nx1, -]
    prob.H=[];              % [Nx1, -]
    
end