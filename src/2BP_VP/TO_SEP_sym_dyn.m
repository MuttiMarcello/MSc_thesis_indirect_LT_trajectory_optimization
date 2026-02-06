function TO_SEP_sym_dyn(prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP variable parameters

% TO symbolic sys dynamics

    y=sym('y',[14,1]);
    assume(y,'real');
    
    rr=y(1:3);
    vv=y(4:6);
    m=y(7);
    llr=y(8:10);
    llv=y(11:13);
    lm=y(14);
    
    U=sym('U',[1,4]); % to be substituted
    assume(U,'real');
    
    LU=U(1);
    TU=U(2);
    MU=U(3);
    cf=U(4);% =1e-6
    
    syms g0; % to be substituted
    assume(g0,'real');
    
    MP=sym('MP',[3,5]);
    assume(MP,'real');
    
    ap=MP(1,:);
    bp=MP(2,:);
    cp=MP(3,:);
    
    %-P=MED-Dyn-computed-inside-F-and-DyF----------------------------------
    
    P=dot(norm(rr).^(0:4),cp);  % [W]
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]

    dP=dot(norm(rr).^(0:3),(1:4).*cp(2:end));
    dT=dot(P.^(0:3),(1:4).*ap(2:end));
    dI=dot(P.^(0:3),(1:4).*bp(2:end));    
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    I=I/TU;               % [-]

    tr=cf*TU^2/(MU*LU)*dT*dP*rr/norm(rr);
    ir=1/TU*dI*dP*rr/norm(rr);
    c=I*TU^2/LU*g0;
    
    % insert S and u dependency in EO/FO
    u=1;
    
    FF_med=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3+u*norm(llv)/m*tr+u*lm/(c)*(tr-T/I*ir);
        -llr;
        -u*norm(llv)*T/m^2];

    A_med=jacobian(FF_med,y);

    %-P=MAX-Dyn-computed-inside-F-and-DyF----------------------------------

    P=sym(prob.Plim(2));
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s] 
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    I=I/TU;               % [-]
    c=I*TU^2/LU*g0;
    
    % insert S and u dependency in EO/FO
    u=1;
    
    FF_max=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        -u*norm(llv)*T/m^2];

    A_max=jacobian(FF_max,y);

    dir_name='TO_dyn';

    if ~isfolder(dir_name)
        mkdir(dir_name)
        addpath(dir_name)
    end

    matlabFunction(FF_med,A_med,'vars',{y,MP,U,g0},'file',fullfile(dir_name,'TO_med_dz'));
    matlabFunction(FF_max,A_max,'vars',{y,MP,U,g0},'file',fullfile(dir_name,'TO_max_dz'));

end