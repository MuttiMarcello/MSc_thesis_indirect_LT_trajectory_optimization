function CR3BP_sym_dyn()
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% TO symbolic sys dynamics

    y=sym('y',[14,1]);
    assume(y,'real');
    
    rr=y(1:3);
    rx=y(1);
    ry=y(2);
    rz=y(3);

    vv=y(4:6);
    vx=vv(1);
    vy=vv(2);

    m=y(7);
    llr=y(8:10);
    llv=y(11:13);
    lm=y(14);

    syms epsilon
    assume(epsilon,'real');
    
    U=sym('U',[1,3]); % to be substituted
    assume(U,'real');
    
    T=U(1);
    c=U(2);
    mu=U(3);
    
    r1=sqrt((rx+mu)^2+ry^2+rz^2);
    r2=sqrt((rx+mu-1)^2+ry^2+rz^2);

    gr=[rx-(1-mu)*(rx+mu)/r1^3-mu*(rx+mu-1)/r2^3;
        ry-(1-mu)*ry/r1^3-mu*ry/r2^3;
        -(1-mu)*rz/r1^3-mu*rz/r2^3];

    hv=[2*vy; -2*vx; 0];

    G11=1-(1-mu)/r1^3+3*(1-mu)*(rx+mu)^2/r1^5-mu/r2^3+3*mu*(rx+mu-1)^2/r2^5;
    G22=1-(1-mu)/r1^3+3*(1-mu)*ry^2/r1^5-mu/r2^3+3*mu*ry^2/r2^5;
    G33=-(1-mu)/r1^3+3*(1-mu)*rz^2/r1^5-mu/r2^3+3*mu*rz^2/r2^5;
    G12=3*(1-mu)*(rx+mu)*ry/r1^5+3*mu*(rx+mu-1)*ry/r2^5;
    G13=3*(1-mu)*(rx+mu)*rz/r1^5+3*mu*(rx+mu-1)*rz/r2^5;
    G23=3*(1-mu)*ry*rz/r1^5+3*mu*ry*rz/r2^5;

    G=[G11, G12, G13; G12, G22, G23; G13, G23, G33];

    H=[0 2 0; -2 0 0; 0 0 0];

    dir_name='CR3BP_Dyn';

    if ~isfolder(dir_name)
        mkdir(dir_name)
        addpath(dir_name)
    end

    %-u_on-Dyn-------------------------------------------------------------
     
    FF=[vv;
        gr+hv-llv/norm(llv)*T/m;
        -T/c;
        -G.'*llv;
        -llr-H.'*llv;
        -norm(llv)*T/m^2];
    
    Aa=sym('Aa',[14,14]);
    Aa(:)=0;

    Aa(1:3,4:6)=eye(3);
    Aa(4:6,1:3)=G;
    Aa(4:6,4:6)=H;
    Aa(4:6,7)=llv/norm(llv)*T/m^2;
    Aa(4:6,11:13)=-T/m*(eye(3)/norm(llv)-llv*llv.'/norm(llv)^3);
    Aa(8:10,1:3)=-jacobian(G.'*llv,rr);
    Aa(8:10,11:13)=-G.';
    Aa(11:13,8:10)=-eye(3);
    Aa(11:13,11:13)=-H.';
    Aa(14,7)=2*norm(llv)*T/m^3;
    Aa(14,11:13)=-llv.'/norm(llv)*T/m^2;

    matlabFunction(FF,Aa,'vars',{y,U},'file',fullfile(dir_name,'u_on_dz'));

    %-u_off-Dyn------------------------------------------------------------
    
    FF=[vv;
        gr+hv;
        0;
        -G.'*llv;
        -llr-H.'*llv;
        0];
    
    Aa=sym('Aa',[14,14]);
    Aa(:)=0;

    Aa(1:3,4:6)=eye(3);
    Aa(4:6,1:3)=G;
    Aa(4:6,4:6)=H;
    Aa(8:10,1:3)=-jacobian(G.'*llv,rr);
    Aa(8:10,11:13)=-G.';
    Aa(11:13,8:10)=-eye(3);
    Aa(11:13,11:13)=-H.';

    matlabFunction(FF,Aa,'vars',{y,U},'file',fullfile(dir_name,'u_off_dz'));

    %-u_med-Dyn------------------------------------------------------------
    
    S=-norm(llv)*c/m-lm+1;
    u=(epsilon-S)/(2*epsilon);

    FF=[vv;
        gr+hv-llv/norm(llv)*u*T/m;
        -u*T/c;
        -G.'*llv;
        -llr-H.'*llv;
        -norm(llv)*u*T/m^2];
    
    Aa=sym('Aa',[14,14]);
    Aa(:)=0;

    Aa(1:3,4:6)=eye(3);
    Aa(4:6,1:3)=G;
    Aa(4:6,4:6)=H;
    Aa(4:6,7)=llv/norm(llv)*u*T/m^2+llv*c*T/(2*epsilon*m^3);
    Aa(4:6,11:13)=-u*T/m*(eye(3)/norm(llv)-llv*llv.'/norm(llv)^3)-llv*llv.'/norm(llv)^2*c*T/(2*epsilon*m^2);
    Aa(4:6,14)=-llv/norm(llv)*T/(2*epsilon*m);
    Aa(7,7)=norm(llv)*T/(2*epsilon*m^2);
    Aa(7,11:13)=-llv.'/norm(llv)*T/(2*epsilon*m);
    Aa(7,14)=-T/(2*epsilon*c);
    Aa(8:10,1:3)=-jacobian(G.'*llv,rr);
    Aa(8:10,11:13)=-G.';
    Aa(11:13,8:10)=-eye(3);
    Aa(11:13,11:13)=-H.';
    Aa(14,7)=2*norm(llv)*u*T/m^3+norm(llv)^2*c*T/(2*epsilon*m^4);
    Aa(14,11:13)=-llv.'/norm(llv)*u*T/m^2-llv.'*c*T/(2*epsilon*m^3);
    Aa(14,14)=-norm(llv)*T/(2*epsilon*m^2);

    matlabFunction(FF,Aa,'vars',{y,U,epsilon},'file',fullfile(dir_name,'u_med_dz'));

end