function [ll0] = ACT(prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters  
    
    x0=[prob.x0; 1];

    rr=x0(1:3);
    vv=x0(4:6);
    m=x0(7);

    v=norm(vv);
    
    a=rad2deg(unifrnd(-10, 10));
    da=rad2deg(unifrnd(-5, 5));
    b=rad2deg(unifrnd(-1, 1));
    db=rad2deg(unifrnd(-0.1, 0.1));

    S=unifrnd(-1, 0)+(prob.isFO-1);
    dS=unifrnd(-0.01, 0.01);
    
    lm=1;

    T=prob.param(1);
    c=prob.param(2);
    mu=prob.param(3);
    
    u=1;

    uup=[cos(a)*cos(b);
        sin(a)*cos(b);
        sin(b)];

    duup=[-sin(a)*da*cos(b)-cos(a)*sin(b)*db;
         cos(a)*da*cos(b)-sin(a)*sin(b)*db;
         cos(b)*db];

    lv=-m/c*(S+lm-prob.isFO);

    dm=-u*T/c;
    dlm=-u*lv*T/m^2;

    dlv=-dm/c*(S+lm)-m/c*(dS+dlm);

    vu=vv/v;

    hh=cross(rr,vv);
    h=norm(hh);
    hu=hh/h;

    R=[vu, cross(hu,vu), hu];
    uu=R*uup;

    llv=-lv*uu;

    rx=rr(1);
    ry=rr(2);
    rz=rr(3);
    vx=vv(1);
    vy=vv(2);

    r1=sqrt((rx+mu)^2+ry^2+rz^2);
    r2=sqrt((rx+mu-1)^2+ry^2+rz^2);

    gr=[rx-(1-mu)*(rx+mu)/r1^3-mu*(rx+mu-1)/r2^3;
        ry-(1-mu)*ry/r1^3-mu*ry/r2^3;
        -(1-mu)*rz/r1^3-mu*rz/r2^3];

    hv=[2*vy; -2*vx; 0];

    H=zeros(3,3);
    H(1,2)=2;
    H(2,1)=-2;

    dvv=gr+hv;
    dv=dot(dvv,vv)/v;
    dvu=dvv/v-dv*vv/v^2;

    dhh=cross(rr,dvv);
    dh=dot(dhh,hh)/h;
    dhu=dhh/h-dh*hh/h^2;

    dR=[dvu, cross(dhu,vu)+cross(hu,dvu), dhu];
    duu=R*duup+dR*uup;

    dllv=-dlv*uu-lv*duu;
    llr=-dllv-H.'*llv;

    ll0=zeros(7,1);
    ll0(1:3)=llr;
    ll0(4:6)=llv;
    ll0(7)=lm;

end