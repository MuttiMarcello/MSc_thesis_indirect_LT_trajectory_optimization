function [df,J,prob] = TO_ZFP(llt,prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% TO zero-finding problem properties given initial costate
% FO_gamma: zfp shooting function
% FO_jacobian: zfp associated jacobian/gradient

    TU=3.75676967e5; % [s]

    t0=prob.t0;

    ll0=llt(1:7);
    tf_ad=llt(8);  % adimensional
    tf=tf_ad*TU+t0;

    prob.tf_ad=tf_ad;
    prob.tf=tf;

    xx0=prob.x0;

    yy0=[xx0; 1; ll0];

    Phi0=eye(length(yy0));
    vPhi0=reshape(Phi0,[length(yy0)^2,1]);

    [~,zz]=TO_ode87(prob,[0 tf_ad],[yy0; vPhi0]);

    zzf=zz(end,:).';

    df=TO_gamma(zzf,prob);
    J=TO_jacobian(zzf,prob);

    prob.y0=yy0;
    prob.gamma=df;
    prob.jac=J;

end

function df = TO_gamma(zf,prob)

    tf_ad=prob.tf_ad;

    rrf=zf(1:3);
    vvf=zf(4:6);

    llrf=zf(8:10);
    llvf=zf(11:13);
    lmf=zf(14);

    ff=TO_CR3BP_cp(tf_ad,zf,prob);
    ffx=ff(1:7);

    Hf=dot([llrf; llvf; lmf],ffx)+1;

    xtf_ad=prob.xf;

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);

    % trensversality condition
    Om=Hf;

    df=zeros(8,1);
    df(1:6)=[rrf; vvf]-[rrtf; vvtf];
    df(7)=lmf;
    df(8)=Om;
end

function J = TO_jacobian(zf,prob)

    tf_ad=prob.tf_ad;

    T=prob.param(1);
    mu=prob.param(3);

    rrf=zf(1:3);
    rx=rrf(1);
    ry=rrf(2);
    rz=rrf(3);
    mf=zf(7);

    r1=sqrt((rx+mu)^2+ry^2+rz^2);
    r2=sqrt((rx+mu-1)^2+ry^2+rz^2);

    G11=1-(1-mu)/r1^3+3*(1-mu)*(rx+mu)^2/r1^5-mu/r2^3+3*mu*(rx+mu-1)^2/r2^5;
    G22=1-(1-mu)/r1^3+3*(1-mu)*ry^2/r1^5-mu/r2^3+3*mu*ry^2/r2^5;
    G33=-(1-mu)/r1^3+3*(1-mu)*rz^2/r1^5-mu/r2^3+3*mu*rz^2/r2^5;
    G12=3*(1-mu)*(rx+mu)*ry/r1^5+3*mu*(rx+mu-1)*ry/r2^5;
    G13=3*(1-mu)*(rx+mu)*rz/r1^5+3*mu*(rx+mu-1)*rz/r2^5;
    G23=3*(1-mu)*ry*rz/r1^5+3*mu*ry*rz/r2^5;

    G=[G11, G12, G13; G12, G22, G23; G13, G23, G33];

    H=[0 2 0; -2 0 0; 0 0 0];

    llvf=zf(11:13);

    lvf=norm(llvf);

    vPhif=zf(15:210);
    Phif=reshape(vPhif,[14,14]);

    ff=TO_CR3BP_cp(tf_ad,zf,prob);
    ffx=ff(1:7);

    dmf=ffx(7);

    dllrf=ff(8:10);
    dllvf=ff(11:13);
    dlmf=ff(14);

    %-Jacobian-------------------------------------------------------------
    J=zeros(8,8);
    
    %-costate-depenndency--------------------------------------------------
    J(1:6,1:7)=Phif(1:6,8:14);      % phi_r,ll phi_v,ll
    J(7,1:7)=Phif(14,8:14);         % phi_lm,ll

    Dffx=zeros([7,7]);

    Dffx(1:3,:)=Phif(4:6,8:14);
    Dffx(4:6,:)=G*Phif(1:3,8:14)+H*Phif(4:6,8:14)+T/(mf*lvf)*(llvf/mf*Phif(7,8:14)+((llvf*llvf.')/(lvf^2)-eye(3))*Phif(11:13,8:14));

    J(8,1:7)=zf(8:14).'*Dffx+ffx.'*Phif(8:14,8:14);

    
    %-time-dependency------------------------------------------------------
    J(1:6,8)=ffx(1:6);   % dpsi/dtf
    J(7,8)=dlmf;

    dffx=zeros(7,1);

    dffx(1:3)=ffx(4:6);
    dffx(4:6)=G*ffx(1:3)+H*ffx(4:6)-T/(mf*lvf)*((eye(3)-(llvf*llvf.')/(lvf^2))*dllvf-dmf/mf*llvf);

    J(8,8)=zf(8:14).'*dffx+ffx.'*[dllrf; dllvf; dlmf];

end

