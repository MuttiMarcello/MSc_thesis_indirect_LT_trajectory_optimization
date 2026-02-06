function [df,J,prob] = FO_ZFP(ll0,prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% EO/FO zero-finding problem properties given initial costate
% FO_gamma: zfp shooting function
% FO_jacobian: zfp associated jacobian/gradient

    tf_ad=prob.tf_ad;

    xx0=prob.x0;

    yy0=[xx0; 1; ll0];

    Phi0=eye(length(yy0));
    vPhi0=reshape(Phi0,[length(yy0)^2,1]);

    [~,zz] =FO_ode87(prob,[0 tf_ad],[yy0; vPhi0]);

    zzf=zz(end,:).';

    df=FO_gamma(zzf,prob);
    J=FO_jacobian(zzf);

    prob.y0=yy0;
    prob.gamma=df;
    prob.jac=J;

end

function df = FO_gamma(zf,prob)

    rrf=zf(1:3);
    vvf=zf(4:6);

    lmf=zf(14);

    xtf_ad=prob.xf;

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);

    df=zeros(7,1);
    df(1:6)=[rrf; vvf]-[rrtf; vvtf];
    df(7)=lmf;

end

function J = FO_jacobian(zf)

    vPhif=zf(15:210);
    Phif=reshape(vPhif,[14,14]);

    %-Jacobian-initialization----------------------------------------------
    J=zeros(7,7);
    
    %-costate-depenndency--------------------------------------------------
    J(1:6,1:7)=Phif(1:6,8:14);      % phi_r,ll phi_v,ll
    J(7,1:7)=Phif(14,8:14);         % phi_lm,ll

end

