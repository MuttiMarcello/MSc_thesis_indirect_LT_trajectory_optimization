function TO_SEP_sym_dyn(prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP constant parameters

% TO symbolic sys dynamics

    y=sym('y',[14,1]);
    assume(y,'real');
    
    rr=y(1:3);
    vv=y(4:6);
    m=y(7);
    llr=y(8:10);
    llv=y(11:13);
    
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
    
    %-P=MAX-Dyn-computed-inside-F-and-DyF----------------------------------

    P=sym(prob.Plim(2));
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s] 
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    I=I/TU;               % [-]
    c=I*TU^2/LU*g0;
    
    % insert S and u dependency in EO/FO
    u=1;
    
    FF=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        -u*norm(llv)*T/m^2];

    A=jacobian(FF,y);

    dir_name='TO_dyn';

    if ~isfolder(dir_name)
        mkdir(dir_name)
        addpath(dir_name)
    end

    matlabFunction(FF,A,'vars',{y,MP,U,g0},'file',fullfile(dir_name,'TO_max_dz'));

end