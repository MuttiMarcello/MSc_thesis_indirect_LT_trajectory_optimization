function FO_SEP_sym_dyn(prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP constant parameters

% EO/FO symbolic sys dynamics

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

    syms epsilon
    assume(epsilon,'real');
 
    %-P=MAX-Dyn-computed-inside-F-and-DyF----------------------------------

    P=sym(prob.Plim(2));
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]   
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    I=I/TU;               % [-]
    c=I*TU^2/LU*g0;

    %-Pmax-uon-------------------------------------------------------------
     
    FF_uon=[vv;
        -rr/norm(rr)^3-T/m*llv/norm(llv);
        -T/c;
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        -norm(llv)*T/m^2];

    A_uon=jacobian(FF_uon,y);

    %-Pmax-uoff------------------------------------------------------------
    
    FF_uoff=[vv;
        -rr/norm(rr)^3;
        0;
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        0];

    A_uoff=jacobian(FF_uoff,y);

    %-Pmax-umed------------------------------------------------------------
    S=-norm(llv)*c/m-lm+1;
    u=(epsilon-S)/(2*epsilon);
        
    FF_umed=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/c;
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        -u*norm(llv)*T/m^2];

    A_umed=jacobian(FF_umed,y);

    
    dir_name='FO_dyn';

    if ~isfolder(dir_name)
        mkdir(dir_name)
        addpath(dir_name)
    end
    
    matlabFunction(FF_uon, A_uon ,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_uon_dz' ));
    matlabFunction(FF_uoff,A_uoff,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_uoff_dz'));
    matlabFunction(FF_umed,A_umed,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_umed_dz'));
    
   
end