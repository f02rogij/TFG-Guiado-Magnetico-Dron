%% ================================================================
%  dron_escalera_autoguiado.m   (3a Aproximacion - AUTOGUIADO PURO)
%  El dron sigue la autopista en ESCALERA EN ESQUINA guiandose SOLO por el
%  campo magnetico (sin trayectoria de referencia).
%   - Sensado: 4 magnetometros (RM3100) en las puntas de los brazos -> gradiente
%     MEDIDO por diferencias finitas; IMU (LSM6DSO) para actitud y velocidades.
%     Cada sensor a su frecuencia, con retencion (ZOH) entre muestras.
%   - Guiado (ejes del cuerpo): rumbo = direccion de minimo |grad B| (la guiñada
%     sigue al hilo); centrado anulando Bperp; avance a v_des con freno en curva;
%     altura PID a z_ref; INS lateral acotada por el campo.
%   - Cuckoo optimiza 13 ganancias en PLANTA LIMPIA; coste = centrado + altura
%     (desviacion real, solo para medir) + suavidad + potencias.
%   - Casos: autoguiado | autoguiado + perturbaciones.
%% ================================================================
clear; close all; tic;
SEED = 42;
[script_dir, prog] = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; prog = 'dron_escalera_autoguiado'; end
outdir = fullfile(script_dir, 'resultados', prog);
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('Salida en: %s\n', outdir);

%% ── Parametros fisicos ──────────────────────────────────────────────
p.m=13.0; p.g=9.81; p.L=0.64; p.Ixx=0.65; p.Iyy=0.65; p.Izz=1.31; p.Ir=4.7e-3;
p.kf=2.07e-4; p.kM=3.76e-6; p.km_el=0.09; p.kdrot=1.5e-5; p.R=0.07;
p.k1e=1.01; p.k2e=0.23; p.k3e=0.74; p.k4e=0.12; p.l_arm=p.L; p.km=p.km_el;
Mact=[p.kf,p.kf,p.kf,p.kf; 0,-p.kf*p.l_arm,0,p.kf*p.l_arm; p.kf*p.l_arm,0,-p.kf*p.l_arm,0; -p.kdrot,p.kdrot,-p.kdrot,p.kdrot];
p.Mact_inv=Mact\eye(4);
p.wPm = 1e-5;   % peso coste mecanico
p.wPe = 1e-4;   % peso coste electrico
p.wRate = 1.5;  % peso velocidad de actitud (rompe ciclo limite)

%% ── Magnetometro (4 sensores) / navegacion ──────────────────────────
% Sensor: PNI RM3100 (magnetoinductivo). Elegido por ser el de MENOR RUIDO de la
%   seleccion (spec critica en un gradiometro: el gradiente=diferencia/baseline
%   amplifica el ruido). Ruido 10-15 nT RMS, hasta 600 Hz, rango +-800 uT
%   (>> campo aqui, decenas de uT con la Tierra incluida).
mag.sigB=13e-9;       % ruido tipico RM3100 (rango 10-15 nT)
mag.Ts=0.01;          % muestreo a 100 Hz (el RM3100 llega a 600 Hz; basta 100 para el lazo)
nav.b_g=2.5e-6; nav.k_obs=0.8;          % deriva INS lateral; ganancia observador (centrado)
nav.b_z=2.5e-6; nav.k_obs_z=0.8;        % deriva INS vertical; ganancia observador de ALTURA (L_z)
nav.b_vf=1.5e-6; nav.k_obs_v=0.5;       % deriva INS avance; ganancia observador de AVANCE (flujo)
nav.lam_a=1e-13;                        % umbral de observabilidad along: separa recta (altura) de travesaño (avance) [TUNEAR]
l_arm=p.L;                              % baseline del gradiometro = brazo
tilt_max=deg2rad(35);
% IMU: ST LSM6DSO (6 ejes). Elegido por el mejor ruido de la gama general + bajo consumo.
%   La deriva del IMU (sesgo sin calibrar) la acota el campo; por eso vale un IMU barato.
imu.ODR=208;                              % frecuencia de medicion usada [Hz] (LSM6DSO hasta 6.6 kHz)
imu.sig_w=(0.0038*pi/180)*sqrt(imu.ODR/2);% ruido vel. angular [rad/s] (giro 0.0038 deg/s/sqrt(Hz))
imu.sig_v=0.01;                           % ruido vel. lineal estimada (INS, repr.) [m/s] (acel 70 ug/sqrt(Hz))
imu.sig_ang=1e-3;                         % ruido actitud estimada (AHRS, repr.) [rad] (~0.06 deg)

%% ── Geometria ESCALERA EN ESQUINA ──────────────────────────────────
Lc=30; w=4; P=8;
r1x=3:P:21; r2y=9:P:Lc;
cond={};
cond{end+1}=struct('nodes',[-2 -w 0; Lc+w -w 0; Lc+w Lc+2 0],'sig',+1);   % riel exterior
cond{end+1}=struct('nodes',[-2 +w 0; Lc-w +w 0; Lc-w Lc+2 0],'sig',-1);   % riel interior
for k=1:numel(r1x), cond{end+1}=struct('nodes',[r1x(k) -w 0; r1x(k) +w 0],'sig',(-1)^k); end %#ok<SAGROW>
for k=1:numel(r2y), cond{end+1}=struct('nodes',[Lc-w r2y(k) 0; Lc+w r2y(k) 0],'sig',(-1)^(k+1)); end %#ok<SAGROW>
wires.cond=cond; wires.I=200; wires.track=[-2 0 0; Lc 0 0; Lc Lc+2 0];

%% ── Consignas de autoguiado ─────────────────────────────────────────
v_des=0.30;    % VELOCIDAD DE AVANCE fija (parametro de control, no optimizable)
z_ref=5.0;     % CRITERIO DE ALTURA (cota constante)
Ltrack=0; for s=1:size(wires.track,1)-1, Ltrack=Ltrack+norm(wires.track(s+1,:)-wires.track(s,:)); end
Tfin=1.30*Ltrack/v_des; Ts=0.005; t_sim=(0:Ts:Tfin)'; N_sim=length(t_sim);
mag_step=max(1,round(mag.Ts/Ts));
imu_step=max(1,round((1/imu.ODR)/Ts));   % retencion IMU (ZOH); =1 si el IMU va al ritmo del lazo
% estado inicial: ya en crucero (v_des sobre el eje x), 0.5 m descentrado para ver el centrado
x0=[-2; v_des; 0.5; 0; z_ref; 0; 0;0;0;0; 0;0];
fprintf('=== Escalera AUTOGUIADO | v_des=%.2f m/s | z_ref=%.1f m | Lpista=%.1f m | Tfin=%.0f s ===\n\n',v_des,z_ref,Ltrack,Tfin);

%% ── Cuckoo: optimizacion sobre PLANTA LIMPIA (13 ganancias) ─────────
%      k_c   kv_f  kv_l  k_fr  kpz  kiz  kdz  kpphi kdphi kpth kdth kppsi kdpsi
lb=[  0.10, 0.20, 0.20, 0.00, 10.0,0.01, 5.0, 40.0, 6.0, 40.0,6.0, 47.0, 5.0];
ub=[  3.00, 5.00, 5.00, 5.00, 70.0,10.0,120.0,260.0,26.0,260.0,26.0,335.0,50.0];
dim=13;
obj=@(v) obj_auto(v,p,t_sim,Ts,N_sim,x0,wires,mag,nav,v_des,z_ref,l_arm,tilt_max,mag_step);
n_pop=30; MaxGen=90; pa=0.25; alpha0=0.1*(ub-lb); lambda=1.5;
rng(SEED);
[best,fval,hi,hf,~]=cuckoo_optimize(obj,lb,ub,dim,n_pop,MaxGen,pa,alpha0,lambda,'autoguiado puro');
G=vec2gains(best); print_gains('GANANCIAS (autoguiado puro)',G,fval);

%% ── Simulaciones finales: 2 casos ───────────────────────────────────
fprintf('\nSimulando casos...\n');
Sa=run_sim(G,p,t_sim,Ts,N_sim,x0,wires,mag,nav,imu,v_des,z_ref,l_arm,tilt_max,mag_step,imu_step,false,true);  % autoguiado
Sp=run_sim(G,p,t_sim,Ts,N_sim,x0,wires,mag,nav,imu,v_des,z_ref,l_arm,tilt_max,mag_step,imu_step,true, true);  % + perturbaciones
rep('Autoguiado         ',Sa,z_ref);
rep('Con perturbaciones ',Sp,z_ref);

%% ── Salidas ─────────────────────────────────────────────────────────
P_=@(n) fullfile(outdir,n);
plot_convergencia(hi,hf,fval,'Cuckoo (autoguiado)',P_('evolucion_cuckoo_autoguiado.png'));
plot_results(Sa,t_sim,z_ref,v_des,Tfin,'Autoguiado magnetico | Escalera',        P_('caso_autoguiado.png'));
plot_results(Sp,t_sim,z_ref,v_des,Tfin,'Autoguiado + perturbaciones | Escalera', P_('caso_perturbaciones_autoguiado.png'));
plot_observabilidad(Sa,t_sim,z_ref,v_des,'Observabilidad - Autoguiado',P_('observabilidad_autoguiado.png'));
plot_3d(Sa,wires,'autoguiado',          P_('trayectoria_3D_autoguiado.png'));
plot_3d(Sp,wires,'con perturbaciones',  P_('trayectoria_3D_perturbaciones_autoguiado.png'));
fprintf('\nCompletado.  Tiempo total: %.1f s\n', toc);


%% ================================================================
%%  FUNCIONES LOCALES
%% ================================================================
function rep(nombre,S,z_ref)
    rec=sum(sqrt(diff(S.X(:,1)).^2+diff(S.X(:,3)).^2+diff(S.X(:,5)).^2));
    vmed=mean(sqrt(S.X(:,2).^2+S.X(:,4).^2));
    fprintf('  %s | recorrido=%.1f m | desv perp RMS(real)=%.3f m | alt RMS err=%.3f m | v media=%.2f m/s | J=%.3f\n', ...
        nombre,rec,sqrt(mean(S.dperp.^2)),sqrt(mean((z_ref-S.X(:,5)).^2)),vmed,S.J(end));
end
function paso=levy_flight(dim,lambda)
    sn=gamma(1+lambda)*sin(pi*lambda/2); sd=gamma((1+lambda)/2)*lambda*2^((lambda-1)/2);
    sigma=(sn/sd)^(1/lambda); u=randn(1,dim)*sigma; v=randn(1,dim); paso=u./(abs(v).^(1/lambda));
end
function [bn,bf,hi,hf,gen]=cuckoo_optimize(obj,lb,ub,dim,n_pop,MaxGen,pa,alpha0,lambda,label)
    alpha=alpha0; nidos=lb+rand(n_pop,dim).*(ub-lb); fitness=zeros(1,n_pop);
    parfor i=1:n_pop, fitness(i)=obj(nidos(i,:)); end %#ok<PFBNS>
    [fitness,idx]=sort(fitness); nidos=nidos(idx,:); bn=nidos(1,:); bf=fitness(1);
    hi=zeros(1,MaxGen); hf=zeros(1,MaxGen); stall=0;
    fprintf('--- Cuckoo [%s] ---  J_ini=%.4f\n',label,bf);
    for gen=1:MaxGen
        for i=1:n_pop
            nu=nidos(i,:)+alpha.*levy_flight(dim,lambda).*randn(1,dim); nu=min(max(nu,lb),ub);
            Qi=obj(nu); j=randi(n_pop); if Qi<fitness(j), nidos(j,:)=nu; fitness(j)=Qi; end
        end
        [fitness,idx]=sort(fitness); nidos=nidos(idx,:);
        n_ab=round(pa*n_pop); idx_ab=(n_pop-n_ab+1):n_pop;
        for k=idx_ab, nidos(k,:)=lb+rand(1,dim).*(ub-lb); end
        Qab=zeros(1,numel(idx_ab)); parfor m=1:numel(idx_ab), Qab(m)=obj(nidos(idx_ab(m),:)); end %#ok<PFBNS>
        fitness(idx_ab)=Qab; [fitness,idx]=sort(fitness); nidos=nidos(idx,:);
        if fitness(1)<bf, bf=fitness(1); bn=nidos(1,:); stall=0; else, stall=stall+1; end
        hi(gen)=gen; hf(gen)=bf;
        if mod(gen,10)==0||gen==1, fprintf('  %-4d  %.4f\n',gen,bf); end
        if stall==15, alpha=alpha*2; end
        if stall>=30, fprintf('  Parada gen %d\n',gen); break; end
    end
    hi=hi(1:gen); hf=hf(1:gen); fprintf('  J_final=%.6f\n\n',bf);
end
function print_gains(t,G,f)
    fprintf('=== %s ===\n',t);
    fprintf('  k_c=%.3f kv_f=%.3f kv_l=%.3f k_freno=%.3f\n',G.k_c,G.kv_f,G.kv_l,G.k_freno);
    fprintf('  kpz=%.3f kiz=%.4f kdz=%.3f\n',G.kpz,G.kiz,G.kdz);
    fprintf('  kpphi=%.2f kdphi=%.2f kptheta=%.2f kdtheta=%.2f kppsi=%.2f kdpsi=%.2f | J=%.4f\n\n',G.kpphi,G.kdphi,G.kptheta,G.kdtheta,G.kppsi,G.kdpsi,f);
end
function G=vec2gains(v)
    G.k_c=v(1);G.kv_f=v(2);G.kv_l=v(3);G.k_freno=v(4);
    G.kpz=v(5);G.kiz=v(6);G.kdz=v(7);
    G.kpphi=v(8);G.kdphi=v(9);G.kptheta=v(10);G.kdtheta=v(11);G.kppsi=v(12);G.kdpsi=v(13);
end
function B=B_segment(r_obs,p1,p2,I)
    mu0=4*pi*1e-7; r_obs=r_obs(:);p1=p1(:);p2=p2(:); lv=p2-p1;L=norm(lv); if L<1e-12,B=zeros(3,1);return;end
    lhat=lv/L; r1=r_obs-p1;r2=r_obs-p2; d_vec=r1-(r1'*lhat)*lhat; d=norm(d_vec); if d<1e-9,B=zeros(3,1);return;end
    cos1=(r1'*lhat)/norm(r1);cos2=(r2'*lhat)/norm(r2); B=(mu0*I)/(4*pi*d)*(cos1-cos2)*cross(lhat,d_vec/d);
end
function B=B_total(r_obs,wires)
    B=zeros(3,1);
    for c=1:numel(wires.cond)
        nd=wires.cond{c}.nodes; Ic=wires.cond{c}.sig*wires.I;
        for s=1:(size(nd,1)-1), B=B+B_segment(r_obs,nd(s,:),nd(s+1,:),Ic); end
    end
end
function y=nz(on,sig,rs)   % ruido gaussiano condicional (0 si no hay ruido)
    if on, y=sig*randn(rs,1,1); else, y=0; end
end
function d=dist_to_track(r,track)
    r=r(:); d=inf;
    for s=1:(size(track,1)-1)
        p1=track(s,:)';p2=track(s+1,:)';lv=p2-p1;L=norm(lv); if L<1e-12,continue;end
        tt=max(0,min(1,dot(r-p1,lv)/L^2)); f=p1+tt*lv;
        dd=hypot(r(1)-f(1),r(2)-f(2));    % distancia HORIZONTAL al eje
        if dd<d, d=dd; end
    end
end
%% ── SENSADO: 4 magnetometros -> gradiente, rumbo, centrado, altura, B0
function m=sense(c,phi,th,psi,l,wires,sigB,rs)
    cphi=cos(phi);sphi=sin(phi);cth=cos(th);sth=sin(th);cps=cos(psi);sps=sin(psi);
    R=[cps*cth, cps*sth*sphi-sps*cphi, cps*sth*cphi+sps*sphi;
       sps*cth, sps*sth*sphi+cps*cphi, sps*sth*cphi-cps*sphi;
       -sth,    cth*sphi,              cth*cphi];
    c=c(:); e1=R*[1;0;0]; e2=R*[0;1;0];
    B1=B_total(c+l*e1,wires); B3=B_total(c-l*e1,wires);
    B2=B_total(c+l*e2,wires); B4=B_total(c-l*e2,wires);
    if sigB>0
        B1=B1+sigB*randn(rs,3,1);B2=B2+sigB*randn(rs,3,1);B3=B3+sigB*randn(rs,3,1);B4=B4+sigB*randn(rs,3,1);
    end
    g1=(B1-B3)/(2*l); g2=(B2-B4)/(2*l);   % dB/d(eje x cuerpo), dB/d(eje y cuerpo) en mundo
    B0=(B1+B2+B3+B4)/4;
    a11=g1'*g1; a12=g1'*g2; a22=g2'*g2;
    % rumbo: direccion (plano horizontal del cuerpo) que minimiza |a*g1+b*g2|
    ang=0.5*atan2(2*a12,a11-a22); wdir=ang+pi/2;
    u=[cos(wdir);sin(wdir)]; if u(1)<0,u=-u;end   % sentido de avance (a>0)
    beta=atan2(u(2),u(1));                  % angulo del hilo en ejes del cuerpo
    % centrado: componente perpendicular (eje y cuerpo) del campo y su sensibilidad
    B0b=R'*B0; g2b=R'*g2;
    Bperp=B0b(2); dBperp=g2b(2);
    lam=(sigB/(2*l))^2+1e-16;
    e_perp=Bperp*dBperp/(dBperp^2+lam); e_perp=max(min(e_perp,5),-5);
    w_obs=dBperp^2/(dBperp^2+lam);
    % altura: gradiente vertical reconstruido por Maxwell (div B=0, rot B=0) -> escala L_z ~ altura
    G1=R'*g1; G2=R'*g2; G3=[G1(3);G2(3);-(G1(1)+G2(2))];   % col z del tensor de gradiente (ejes cuerpo)
    dnorm_dz=(B0b'*G3)/max(norm(B0b),1e-15);                % d|B|/dz (cuerpo-z ~ vertical en vuelo nivelado)
    m.Lz=norm(B0b)/max(abs(dnorm_dz),1e-12);                % L_z = |B|/|d|B|/dz| ~ altura sobre los hilos
    m.beta=beta; m.e_perp=e_perp; m.w_obs=w_obs; m.dBperp=dBperp;
    m.B0=B0; m.g1=g1; m.a11=a11;                            % para dB/dt (avance) y peso de observabilidad
end
function [om,Im,Vm,Pm,Pe]=actuator_model(F,t2,t3,t4,xd_,yd_,zd_,phi,th,psi,dphi,dth,dps,p)
    w2=max(p.Mact_inv*[F;t2;t3;t4],0);om=sqrt(w2);Im=p.kdrot*w2/p.km;Vm=p.R*Im+p.km*om;
    Rez=[cos(psi)*sin(th)+sin(psi)*sin(phi)*cos(th);sin(psi)*sin(th)-cos(psi)*sin(phi)*cos(th);cos(phi)*cos(th)];
    Pm=F*(Rez'*[xd_;yd_;zd_])+t2*dphi+t3*dth+t4*dps;Pe=p.R*sum(Im.^2);
end
function dx=sys_deriv(x,F,t2,t3,t4,mp,Ix,Iy,Iz,p,pert,vf)
    xdot=x(2);ydot=x(4);zdot=x(6);phi=x(7);phidot=x(8);th=x(9);thdot=x(10);psi=x(11);psidot=x(12);
    cp=cos(psi);sp=sin(psi);
    xdd=p.g*(th*cp+phi*sp);ydd=p.g*(th*sp-phi*cp);zdd=F/mp-p.g;phidd=t2/Ix;thdd=t3/Iy;psidd=t4/Iz;
    if pert
        ve=[xdot;ydot;zdot]+0.05*[phidot;thdot;0]-vf;
        Fa1=-p.k1e*sqrt(ve(1)^2+ve(2)^2)*ve(1)-p.k2e*abs(ve(1))*ve(1);
        Fa2=-p.k1e*sqrt(ve(1)^2+ve(2)^2)*ve(2)-p.k2e*abs(ve(2))*ve(2);
        Fa3=-p.k3e*(ve(1)^2+ve(2)^2)-p.k4e*abs(ve(3))*ve(3);
        Fae=[Fa1;Fa2;Fa3]-0.28*sqrt(abs(F))*vf;Mae=-0.05*cross([0;0;1],Fae);
        xdd=xdd+Fae(1)/mp;ydd=ydd+Fae(2)/mp;zdd=zdd+Fae(3)/mp;phidd=phidd+Mae(1)/Ix;thdd=thdd+Mae(2)/Iy;psidd=psidd+Mae(3)/Iz;
    end
    dx=[xdot;xdd;ydot;ydd;zdot;zdd;phidot;phidd;thdot;thdd;psidot;psidd];
end
function dx=fast_deriv(x,F,t2,t3,t4,g_,m_,Ix_,Iy_,Iz_)
    cp=cos(x(11));sp=sin(x(11));
    dx=[x(2); g_*(x(9)*cp+x(7)*sp); x(4); g_*(x(9)*sp-x(7)*cp); x(6); F/m_-g_; x(8); t2/Ix_; x(10); t3/Iy_; x(12); t4/Iz_];
end
%% ── run_sim_fast: planta LIMPIA (sin ruido, sin perturbaciones) -> J ─
function [J,div]=run_sim_fast(G,p,t,Ts,N,x0,wires,~,nav,v_des,z_ref,l,tilt_max,mag_step)
    g_=p.g;m_=p.m;Ix_=p.Ixx;Iy_=p.Iyy;Iz_=p.Izz;
    Ts2=Ts*.5;Ts6=Ts/6;ilim=50;a_b=0.10;a_r=0.10;
    x=x0; int_ez=0; off_hat=0; beta_f=0; beta_prev=0; rate_f=0; Jacc=0; div=false; J=0;
    e_perp=0; w_obs=0; ctr=0;
    for k=1:N
        tk=t(k); phi=x(7);th=x(9);psi=x(11);
        if ctr==0
            mm=sense([x(1);x(3);x(5)],phi,th,psi,l,wires,0,[]);
            beta=mm.beta; e_perp=mm.e_perp; w_obs=mm.w_obs;
            beta_f=beta_f+a_b*(beta-beta_f);
            rate=(beta_f-beta_prev)/(mag_step*Ts); beta_prev=beta_f; rate_f=rate_f+a_r*(rate-rate_f);
        end
        ctr=mod(ctr+1,mag_step);
        vx=x(2);vy=x(4);vz=x(6);phidot=x(8);thdot=x(10);psidot=x(12);
        vfwd=vx*cos(psi)+vy*sin(psi); vlat=-vx*sin(psi)+vy*cos(psi);
        % INS lateral acotada por campo (observador complementario)
        bias=nav.b_g*g_*tk;
        off_hat=off_hat+Ts*(vlat+bias)-Ts*nav.k_obs*w_obs*(off_hat-e_perp);
        % consignas
        v_eff=v_des/(1+G.k_freno*abs(rate_f));
        thd=min(max(G.kv_f*(v_eff-vfwd),-tilt_max),tilt_max);
        phd=min(max(-G.kv_l*(-G.k_c*off_hat-vlat),-tilt_max),tilt_max);
        e_z=z_ref-x(5); int_ez=max(min(int_ez+e_z*Ts,ilim),-ilim);
        F=max(m_*g_+G.kpz*e_z+G.kiz*int_ez-G.kdz*vz,0.05);
        t2=G.kpphi*(phd-phi)-G.kdphi*phidot; t3=G.kptheta*(thd-th)-G.kdtheta*thdot; t4=G.kppsi*beta_f-G.kdpsi*psidot;
        dperp=dist_to_track([x(1);x(3);x(5)],wires.track);
        if dperp>10||abs(phi)>1.2||abs(th)>1.2||~isfinite(x(1)),div=true;return;end
        w2=max(p.Mact_inv*[F;t2;t3;t4],0);Im2=p.kdrot*w2/p.km;
        Rez=[cos(psi)*sin(th)+sin(psi)*sin(phi)*cos(th);sin(psi)*sin(th)-cos(psi)*sin(phi)*cos(th);cos(phi)*cos(th)];
        Pm=F*(Rez'*[vx;vy;vz])+t2*phidot+t3*thdot+t4*psidot;Pe=p.R*sum(Im2.^2);
        Jacc=Jacc+Ts*(tk*(dperp^2+e_z^2)+p.wPm*Pm^2+p.wPe*Pe+p.wRate*(phidot^2+thdot^2));
        f1=fast_deriv(x,F,t2,t3,t4,g_,m_,Ix_,Iy_,Iz_);
        f2=fast_deriv(x+Ts2*f1,F,t2,t3,t4,g_,m_,Ix_,Iy_,Iz_);
        f3=fast_deriv(x+Ts2*f2,F,t2,t3,t4,g_,m_,Ix_,Iy_,Iz_);
        f4=fast_deriv(x+Ts*f3,F,t2,t3,t4,g_,m_,Ix_,Iy_,Iz_);
        x=x+Ts6*(f1+2*f2+2*f3+f4);
    end
    J=Jacc;
end
%% ── run_sim: simulacion completa (con ruido / perturbaciones), guarda todo
function S=run_sim(G,p,t,Ts,N,x0,wires,mag,nav,imu,v_des,z_ref,l,tilt_max,mag_step,imu_step,with_perturb,use_noise)
    if with_perturb, mp=1.10*p.m;Ixp=1.10*p.Ixx;Iyp=1.10*p.Iyy;Izp=1.10*p.Izz;vf=[0.25;0.15;0.0];
    else, mp=p.m;Ixp=p.Ixx;Iyp=p.Iyy;Izp=p.Izz;vf=[0;0;0]; end
    if use_noise, rs=RandStream('mt19937ar','Seed',7); sigB=mag.sigB; else, rs=[]; sigB=0; end
    ilim=50;a_b=0.10;a_r=0.10; x=x0; int_ez=0; off_hat=0; beta_f=0; beta_prev=0; rate_f=0; Jacc=0;
    e_perp=0; w_obs=0; dBperp=0; ctr=0; imu_ctr=0;
    % observadores altura/avance (correccion magnetica) + calibracion L_z
    z_hat=x0(5); b_fhat=0; Lz_cal=0; cal_done=false; B0_prev=zeros(3,1); dBdt=zeros(3,1);
    mg1=zeros(3,1); ma11=0; mLz=x0(5);   % retenidos entre muestras del magnetometro
    % medidas IMU retenidas (init = estado inicial)
    phi_m=x0(7);th_m=x0(9);psi_m=x0(11);p_m=x0(8);q_m=x0(10);r_m=x0(12);vx_m=x0(2);vy_m=x0(4);vz_m=x0(6);
    X=zeros(N,12);U=zeros(N,4);Ang=zeros(N,2);Jo=zeros(N,1);Pmo=zeros(N,1);Peo=zeros(N,1);
    betao=zeros(N,1);offo=zeros(N,1);epo=zeros(N,1);dpo=zeros(N,1);vfo=zeros(N,1);veo=zeros(N,1);wo=zeros(N,1);dBo=zeros(N,1);
    zho=zeros(N,1);zmo=zeros(N,1);vfho=zeros(N,1);vfmo=zeros(N,1);walo=zeros(N,1);
    for k=1:N
        tk=t(k); phi=x(7);th=x(9);psi=x(11);
        if ctr==0
            m=sense([x(1);x(3);x(5)],phi,th,psi,l,wires,sigB,rs);
            beta_f=beta_f+a_b*(m.beta-beta_f);
            rate=(beta_f-beta_prev)/(mag_step*Ts); beta_prev=beta_f; rate_f=rate_f+a_r*(rate-rate_f);
            e_perp=m.e_perp; w_obs=m.w_obs; dBperp=m.dBperp; mg1=m.g1; ma11=m.a11; mLz=m.Lz;
            if cal_done, dBdt=(m.B0-B0_prev)/(mag_step*Ts); else, Lz_cal=m.Lz; cal_done=true; end
            B0_prev=m.B0;
        end
        ctr=mod(ctr+1,mag_step);
        % --- IMU (LSM6DSO): muestreo a su ODR con RETENCION (ZOH) + ruido (solo validacion) ---
        if imu_ctr==0
            phi_m=phi+nz(use_noise,imu.sig_ang,rs); th_m=th+nz(use_noise,imu.sig_ang,rs); psi_m=psi+nz(use_noise,imu.sig_ang,rs);
            p_m=x(8)+nz(use_noise,imu.sig_w,rs); q_m=x(10)+nz(use_noise,imu.sig_w,rs); r_m=x(12)+nz(use_noise,imu.sig_w,rs);
            vx_m=x(2)+nz(use_noise,imu.sig_v,rs); vy_m=x(4)+nz(use_noise,imu.sig_v,rs); vz_m=x(6)+nz(use_noise,imu.sig_v,rs);
        end
        imu_ctr=mod(imu_ctr+1,imu_step);
        % velocidades en ejes del cuerpo segun el IMU (estimadas)
        vfwd=vx_m*cos(psi_m)+vy_m*sin(psi_m); vlat=-vx_m*sin(psi_m)+vy_m*cos(psi_m);
        bias=nav.b_g*p.g*tk;
        off_hat=off_hat+Ts*(vlat+bias)-Ts*nav.k_obs*w_obs*(off_hat-e_perp);
        % observabilidad complementaria (ma11=|grad B| a lo largo del hilo): avance EN travesaños, altura ENTRE ellos
        w_along=ma11/(ma11+nav.lam_a); w_zalt=nav.lam_a/(ma11+nav.lam_a);
        % ALTURA: L_z (entre travesaños) acota la deriva vertical de la INS
        vz_imu=vz_m+nav.b_z*p.g*tk;
        z_mag=z_ref+(mLz-Lz_cal);
        z_hat=z_hat+Ts*vz_imu-Ts*nav.k_obs_z*w_zalt*(z_hat-z_mag);
        % AVANCE: flujo magnetico dB/dt=G*v (en travesaños) acota la deriva de v_fwd
        v_fwd_mag=(dBdt'*mg1)/(ma11+nav.lam_a);
        v_fwd_imu=vfwd+nav.b_vf*p.g*tk;
        b_fhat=b_fhat+Ts*nav.k_obs_v*w_along*((v_fwd_imu-v_fwd_mag)-b_fhat);
        v_fwd_hat=v_fwd_imu-b_fhat;
        v_eff=v_des/(1+G.k_freno*abs(rate_f));
        thd=min(max(G.kv_f*(v_eff-v_fwd_hat),-tilt_max),tilt_max);
        phd=min(max(-G.kv_l*(-G.k_c*off_hat-vlat),-tilt_max),tilt_max);
        e_z=z_ref-z_hat; int_ez=max(min(int_ez+e_z*Ts,ilim),-ilim);
        F=max(p.m*p.g+G.kpz*e_z+G.kiz*int_ez-G.kdz*vz_imu,0.05);
        t2=G.kpphi*(phd-phi_m)-G.kdphi*p_m; t3=G.kptheta*(thd-th_m)-G.kdtheta*q_m; t4=G.kppsi*beta_f-G.kdpsi*r_m;
        [~,~,~,Pm,Pe]=actuator_model(F,t2,t3,t4,x(2),x(4),x(6),phi,th,psi,x(8),x(10),x(12),p);  % potencias con estado REAL
        dperp=dist_to_track([x(1);x(3);x(5)],wires.track);
        e_z_real=z_ref-x(5);                                                                     % coste con error REAL
        Jacc=Jacc+Ts*(tk*(dperp^2+e_z_real^2)+p.wPm*Pm^2+p.wPe*Pe+p.wRate*(x(8)^2+x(10)^2));
        U(k,:)=[F,t2,t3,t4];Ang(k,:)=[phd,thd];Pmo(k)=Pm;Peo(k)=Pe;Jo(k)=Jacc;
        vfwd_real=x(2)*cos(x(11))+x(4)*sin(x(11));   % velocidad de avance REAL (para grafica)
        betao(k)=beta_f;offo(k)=off_hat;epo(k)=e_perp;dpo(k)=dperp;vfo(k)=vfwd_real;veo(k)=v_eff;wo(k)=w_obs;dBo(k)=dBperp;
        zho(k)=z_hat;zmo(k)=z_mag;vfho(k)=v_fwd_hat;vfmo(k)=v_fwd_mag;walo(k)=w_along;
        f1=sys_deriv(x,F,t2,t3,t4,mp,Ixp,Iyp,Izp,p,with_perturb,vf);
        f2=sys_deriv(x+Ts/2*f1,F,t2,t3,t4,mp,Ixp,Iyp,Izp,p,with_perturb,vf);
        f3=sys_deriv(x+Ts/2*f2,F,t2,t3,t4,mp,Ixp,Iyp,Izp,p,with_perturb,vf);
        f4=sys_deriv(x+Ts*f3,F,t2,t3,t4,mp,Ixp,Iyp,Izp,p,with_perturb,vf);
        x=x+(Ts/6)*(f1+2*f2+2*f3+f4); X(k,:)=x';
    end
    S.X=X;S.U=U;S.Ang=Ang;S.J=Jo;S.Pm=Pmo;S.Pe=Peo;
    S.beta=betao;S.off_hat=offo;S.e_perp=epo;S.dperp=dpo;S.vfwd=vfo;S.v_eff=veo;S.w_obs=wo;S.dBperp=dBo;
    S.z_hat=zho;S.z_mag=zmo;S.v_fwd_hat=vfho;S.v_fwd_mag=vfmo;S.w_along=walo;
end
function J=obj_auto(v,p,t,Ts,N,x0,wires,mag,nav,v_des,z_ref,l,tilt_max,mag_step)
    G=vec2gains(v);[J,div]=run_sim_fast(G,p,t,Ts,N,x0,wires,mag,nav,v_des,z_ref,l,tilt_max,mag_step);
    if div||~isfinite(J),J=1e8;end
end
function plot_convergencia(hi,hf,fval,metodo,fp)
    if numel(hf)<2,return;end;[vm,im]=min(hf);
    fig=figure('Position',[100 80 1000 460],'Color','w');
    semilogy(hi,hf,'b-o','MarkerSize',3,'LineWidth',1.8);hold on;semilogy(im,vm,'r*','MarkerSize',14,'LineWidth',2.5);
    xlabel('Iteracion');ylabel('J (log)');title(['Convergencia - ' metodo]);grid on;
    annotation('textbox',[0.6 0.6 0.28 0.2],'String',sprintf('J_{final}=%.4f',fval),'FitBoxToText','on','BackgroundColor',[1 1 .88],'FontName','Courier New','FontSize',9,'Interpreter','none');
    drawnow;try exportgraphics(fig,fp,'Resolution',300);catch,print(fig,fp,'-dpng','-r300');end;close(fig);
end
function plot_results(S,t,z_ref,v_des,Tfin,titulo,fp)
    x=S.X(:,1);y=S.X(:,3);z=S.X(:,5);phi=S.X(:,7);theta=S.X(:,9);psi=S.X(:,11);
    phd=S.Ang(:,1);thd=S.Ang(:,2);xl=[0 Tfin];
    fig=figure('Position',[20 20 1600 1100],'Color','w','Name',titulo);
    subplot(4,3,1);plot(t,x,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('x[m]');title('Posicion X');
    subplot(4,3,2);plot(t,y,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('y[m]');title('Posicion Y');
    subplot(4,3,3);plot(t,z,'b',t,z_ref*ones(size(t)),'r--','LineWidth',1.2);grid on;xlim(xl);ylabel('z[m]');title('Altura (z_{ref} en rojo)');
    subplot(4,3,4);plot(t,thd*180/pi,'r--',t,theta*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('theta[deg]');title('Cabeceo (avance)');
    subplot(4,3,5);plot(t,phd*180/pi,'r--',t,phi*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('phi[deg]');title('Alabeo (centrado)');
    subplot(4,3,6);plot(t,psi*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('psi[deg]');title('Guinada (sigue al hilo)');
    subplot(4,3,7);plot(t,S.dperp,'b','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('d_{perp}[m]');title('Desviacion perp REAL');
    subplot(4,3,8);plot(t,z_ref-z,'g','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('e_z[m]');title('Error de altura');
    subplot(4,3,9);plot(t,S.vfwd,'b',t,v_des*ones(size(t)),'r--',t,S.v_eff,'m:','LineWidth',1.2);grid on;xlim(xl);ylabel('v[m/s]');title('Avance: real / v_{des} / v_{eff}');
    subplot(4,3,10);plot(t,S.beta*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('beta[deg]');title('Desalineamiento con el hilo');
    subplot(4,3,11);plot(t,S.U(:,2),'b',t,S.U(:,3),'r',t,S.U(:,4),'g','LineWidth',1.0);grid on;xlim(xl);ylabel('[Nm]');xlabel('t[s]');legend('tau2','tau3','tau4');title('Torques');
    subplot(4,3,12);yyaxis left;plot(t,S.U(:,1),'b','LineWidth',1.2);ylabel('F[N]');yyaxis right;plot(t,S.J,'m','LineWidth',1.4);ylabel('J');grid on;xlim(xl);xlabel('t[s]');title(sprintf('Empuje y J (J_{fin}=%.2f)',S.J(end)));
    set(findall(fig,'-property','FontSize'),'FontSize',8);sgtitle(titulo,'FontSize',11,'FontWeight','bold');drawnow;
    try exportgraphics(fig,fp,'Resolution',300);catch,print(fig,fp,'-dpng','-r300');end;close(fig);
end
function plot_observabilidad(S,t,z_ref,v_des,titulo,fp)
    xl=[0 t(end)];fig=figure('Position',[40 30 1100 1120],'Color','w','Name',titulo);
    subplot(4,1,1);plot(t,S.dperp,'b',t,S.off_hat,'r--',t,S.e_perp,'g:','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');
    ylabel('[m]');title('Centrado: real / off\_hat (observador) / e\_perp (magnetico)');legend('real','off\_hat','e\_perp');
    subplot(4,1,2);plot(t,S.X(:,5),'b',t,S.z_hat,'r--',t,S.z_mag,'g:',t,z_ref*ones(size(t)),'k-.','LineWidth',1.2);grid on;xlim(xl);
    ylabel('z [m]');title('Altura: z real / z\_hat (observador) / z\_mag (L_z) / z\_ref');legend('z real','z\_hat','z\_mag','z\_ref');
    subplot(4,1,3);plot(t,S.vfwd,'b',t,S.v_fwd_hat,'r--',t,S.v_fwd_mag,'g:',t,v_des*ones(size(t)),'k-.','LineWidth',1.2);grid on;xlim(xl);ylim([-0.2 0.6]);
    ylabel('v_{avance} [m/s]');title('Avance: real / v\_fwd\_hat (observador) / v\_fwd\_mag (flujo) / v\_des');legend('real','hat','mag','v_{des}');
    subplot(4,1,4);plot(t,S.w_obs,'b',t,S.w_along,'r',t,1-S.w_along,'g','LineWidth',1.2);grid on;xlim(xl);
    ylabel('peso');xlabel('t[s]');title('Observabilidad: centrado / avance (en travesaños) / altura (entre travesaños)');legend('w_{centrado}','w_{avance}','w_{altura}');
    sgtitle(titulo,'FontSize',12,'FontWeight','bold');drawnow;
    try exportgraphics(fig,fp,'Resolution',300);catch,print(fig,fp,'-dpng','-r300');end;close(fig);
end
function plot_3d(S,wires,etiqueta,fp)
    fig=figure('Position',[50 50 760 580],'Color','w');
    ax=axes(fig);hold(ax,'on');grid(ax,'on');view(ax,40,28);
    for c=1:numel(wires.cond)
        nd=wires.cond{c}.nodes;
        plot3(ax,nd(:,1),nd(:,2),nd(:,3),'-','Color',[.6 .6 .6],'LineWidth',1.0,'HandleVisibility','off');
    end
    tr=wires.track; plot3(ax,tr(:,1),tr(:,2),tr(:,3),'k--','LineWidth',1.5);
    plot3(ax,S.X(:,1),S.X(:,3),S.X(:,5),'-','Color',[0 .45 .74],'LineWidth',1.8);
    plot3(ax,S.X(1,1),S.X(1,3),S.X(1,5),'ko','MarkerFaceColor','w','MarkerSize',8);
    plot3(ax,S.X(end,1),S.X(end,3),S.X(end,5),'k^','MarkerFaceColor',[0 .45 .74],'MarkerSize',8);
    allp=[S.X(:,[1 3 5]); tr]; for c=1:numel(wires.cond), allp=[allp; wires.cond{c}.nodes]; end %#ok<AGROW>
    mn=min(allp,[],1); mx=max(allp,[],1); mrg=max(0.05*max(mx-mn),0.6);
    daspect(ax,[1 1 1]); xlim(ax,[mn(1)-mrg mx(1)+mrg]); ylim(ax,[mn(2)-mrg mx(2)+mrg]); zlim(ax,[mn(3)-mrg mx(3)+mrg]);
    xlabel(ax,'x[m]');ylabel(ax,'y[m]');zlabel(ax,'z[m]');
    title(sprintf('Trayectoria 3D autoguiado - %s',etiqueta),'FontWeight','bold');drawnow;
    try exportgraphics(fig,fp,'Resolution',300);catch,print(fig,fp,'-dpng','-r300');end;close(fig);
end
