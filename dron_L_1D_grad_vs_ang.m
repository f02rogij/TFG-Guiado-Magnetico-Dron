%% ================================================================
%  dron_L_1D_grad_vs_ang.m   (3a Aproximacion)
%  ----------------------------------------------------------------
%  OBJETIVO: comparar los DOS metodos de correccion magnetica de la
%  desviacion PERPENDICULAR al movimiento (cross-track) sobre una
%  autopista magnetica en L:
%     - GRADIENTE vectorial de B   (delta = -(Jp.e_B)/(||Jp||^2+lam))
%     - ANGULO del vector B        (delta = ga*r_theta/(ga^2+lam_a))
%  Solo correccion PERP (1D). Sin perturbaciones aero.
%  CASOS:  control (sin correccion)  |  grad  |  ang
%  Coste J sobre el error REAL.  Hipotesis: el metodo por GRADIENTES
%  es mejor (usa magnitud + direccion; el de angulo solo direccion).
%
%  NOTA nombre: MATLAB no ejecuta scripts con guiones en el nombre
%  ('grad-vs-ang' = resta), por eso se usan guiones bajos.
%% ================================================================
clear; close all; tic;
SEED = 42;
[script_dir, prog] = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; prog = 'dron_L_1D_grad_vs_ang'; end
outdir = fullfile(script_dir, 'resultados', prog);   % resultados\<nombre_programa>
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('Salida en: %s\n', outdir);

%% ── Parametros fisicos ──────────────────────────────────────────────
p.m     = 13.0;   p.g     = 9.81;   p.L     = 0.64;
p.Ixx   = 0.65;   p.Iyy   = 0.65;   p.Izz   = 1.31;
p.Ir    = 4.7e-3;
p.kf    = 2.07e-4; p.kM    = 3.76e-6; p.km_el = 0.09;
p.kdrot = 1.5e-5;  p.R     = 0.07;
p.k1e   = 1.01;    p.k2e   = 0.23;
p.k3e   = 0.74;    p.k4e   = 0.12;
p.l_arm = p.L;     p.km    = p.km_el;
Mact = [ p.kf,          p.kf,          p.kf,          p.kf;
         0,            -p.kf*p.l_arm,  0,             p.kf*p.l_arm;
         p.kf*p.l_arm,  0,            -p.kf*p.l_arm,  0;
        -p.kdrot,       p.kdrot,      -p.kdrot,       p.kdrot ];
p.Mact_inv = Mact \ eye(4);
p.wPm = 1e-5;   % peso coste mecanico  (calibrado: termino ~orden de unidad, no perturba)
p.wPe = 1e-4;   % peso coste electrico (calibrado: termino ~orden de unidad, no perturba)
p.wRate = 1.5;  % peso VELOCIDAD de actitud (phidot^2+thetadot^2) -> rompe el ciclo limite (~10% de J)

%% ── Magnetometro / navegacion ───────────────────────────────────────
mag.sigB    = 15e-9;   mag.Ts = 1/600;   mag.dr = 0.05;   mag.sigma_s = 0.5;
nav.b_g     = 2.5e-6;  nav.k_obs = 0.8;

%% ── Cable en L + trayectoria con offset lateral ─────────────────────
wires.nodes = [0,0,0; 15,0,0; 15,15,0];
wires.I     = 200;
off = 3.5;
wp = [0,      -off, 5, 0;
      15+off, -off, 5, 0;
      15+off, -off, 5, 0;
      15+off,  15,  5, 0];
v_avg = 0.23;  t_dwell = 5;
n_seg = size(wp,1)-1;  d_seg = zeros(n_seg,1); T_seg = zeros(n_seg,1);
for s = 1:n_seg
    d_seg(s) = norm(wp(s+1,1:3)-wp(s,1:3));
    if d_seg(s) > 1e-9, T_seg(s) = d_seg(s)/v_avg; else, T_seg(s) = t_dwell; end
end
t_wp = [0; cumsum(T_seg)];  Tfin = t_wp(end);
fprintf('=== L offset %.1f m | v=%.3f m/s | Tfin=%.1f s ===\n\n', off, v_avg, Tfin);

Ts = 0.005;  t_sim = (0:Ts:Tfin)';  N_sim = length(t_sim);
rd_full = get_trajectory(t_sim, wp, t_wp, n_seg);
x0_drone    = [wp(1,1);0;wp(1,2);0;wp(1,3);0;0;0;0;0;wp(1,4);0];
err_max_opt = 30;

fprintf('Precalculando tabla magnetica...\n');
precomp = make_precomp(t_sim, rd_full, mag, wires, Ts);
fprintf('  Listo.\n\n');

%% ── Espacio de busqueda Cuckoo ──────────────────────────────────────
lb = [  0.1,0.1,  1.0,  0.1,0.001,  1.0, 10.0,0.01,  5.0,  40.0, 6.0,  40.0, 6.0,  47.0, 5.0];
ub = [ 15.0,5.0,180.0, 15.0, 20.0,180.0, 70.0,10.0,120.0, 260.0,26.0, 260.0,26.0, 335.0,50.0];   % kpx,kpy a 15 (antes 150): evita el bang-bang de actitud
dim = 15;
% Optimizacion sobre PLANTA LIMPIA (perp sin deriva). Mismas ganancias
% para los tres casos -> comparacion del metodo SIN confundir con ganancias.
obj = @(v) obj_quad(v, p, rd_full, t_sim, Ts, N_sim, x0_drone, err_max_opt, precomp, mag, wires, nav, 'none');
n_pop = 40; MaxGen = 110; pa = 0.25; alpha0 = 0.1*(ub-lb); lambda = 1.5;

rng(SEED);
[best, fval, hi, hf, ~] = cuckoo_optimize(obj, lb, ub, dim, n_pop, MaxGen, pa, alpha0, lambda, 'planta limpia (comun)');
G = vec2gains(best);
print_gains('GANANCIAS comunes (planta limpia)', G, fval);

%% ── Simulaciones finales: 3 casos (mismas ganancias) ────────────────
fprintf('\nSimulando casos...\n');
S_ctrl = run_sim(G, p, t_sim, rd_full, Ts, N_sim, x0_drone, precomp, mag, wires, nav, 'none');
S_grad = run_sim(G, p, t_sim, rd_full, Ts, N_sim, x0_drone, precomp, mag, wires, nav, 'grad');
S_ang  = run_sim(G, p, t_sim, rd_full, Ts, N_sim, x0_drone, precomp, mag, wires, nav, 'ang');
fprintf('  Control (sin correccion): J=%.4f  | desv perp RMS=%.4f m\n', S_ctrl.J(end), rms_(S_ctrl.dp_true));
fprintf('  GRAD   : J=%.4f  | desv perp RMS=%.4f m\n', S_grad.J(end), rms_(S_grad.dp_true));
fprintf('  ANG    : J=%.4f  | desv perp RMS=%.4f m\n', S_ang.J(end),  rms_(S_ang.dp_true));
Jtot=S_grad.J(end); Jmec=p.wPm*Ts*sum(S_grad.Pm.^2); Jele=p.wPe*Ts*sum(S_grad.Pe); Jrate=p.wRate*Ts*sum(S_grad.X(:,8).^2+S_grad.X(:,10).^2);
fprintf('  Desglose J(grad): seguim=%.0f(%.1f%%)  vel_act=%.0f(%.1f%%)  mecanico=%.2f(%.2f%%)  electrico=%.2f(%.2f%%)\n', ...
    Jtot-Jmec-Jele-Jrate,100*(Jtot-Jmec-Jele-Jrate)/Jtot, Jrate,100*Jrate/Jtot, Jmec,100*Jmec/Jtot, Jele,100*Jele/Jtot);

%% ── Salidas ─────────────────────────────────────────────────────────
P = @(n) fullfile(outdir, n);
plot_convergencia(hi, hf, fval, G, 'Cuckoo (planta limpia)', P('evolucion_cuckoo_1D.png'));
plot_results(S_ctrl, t_sim, rd_full, Tfin, 'Control (sin correccion) | L 1D perp', P('caso_control_1D.png'));
plot_results(S_grad, t_sim, rd_full, Tfin, 'Metodo GRADIENTE | L 1D perp',          P('caso_grad_1D.png'));
plot_results(S_ang,  t_sim, rd_full, Tfin, 'Metodo ANGULO | L 1D perp',             P('caso_ang_1D.png'));
plot_compara_1D(S_ctrl, S_grad, S_ang, t_sim, P('comparativa_1D_grad_vs_ang.png'));

fprintf('\nCompletado.  Tiempo total: %.1f s\n', toc);


%% ================================================================
%%  FUNCIONES LOCALES
%% ================================================================
function paso = levy_flight(dim, lambda)
    sn=gamma(1+lambda)*sin(pi*lambda/2); sd=gamma((1+lambda)/2)*lambda*2^((lambda-1)/2);
    sigma=(sn/sd)^(1/lambda); u=randn(1,dim)*sigma; v=randn(1,dim); paso=u./(abs(v).^(1/lambda));
end

function [best_nido, best_fitness, hist_iter, hist_fval, gen] = ...
        cuckoo_optimize(obj, lb, ub, dim, n_pop, MaxGen, pa, alpha0, lambda, label)
    alpha = alpha0;
    nidos = lb + rand(n_pop,dim).*(ub-lb);
    fitness = zeros(1,n_pop);
    parfor i = 1:n_pop, fitness(i) = obj(nidos(i,:)); end %#ok<PFBNS>
    [fitness,idx] = sort(fitness); nidos = nidos(idx,:);
    best_nido = nidos(1,:); best_fitness = fitness(1);
    hist_iter = zeros(1,MaxGen); hist_fval = zeros(1,MaxGen); stall = 0;
    fprintf('--- Cuckoo [%s] ---  J_ini = %.4f\n', label, best_fitness);
    for gen = 1:MaxGen
        for i = 1:n_pop
            nuevo = nidos(i,:) + alpha.*levy_flight(dim,lambda).*randn(1,dim);
            nuevo = min(max(nuevo,lb),ub);
            Qi = obj(nuevo);  j = randi(n_pop);
            if Qi < fitness(j), nidos(j,:)=nuevo; fitness(j)=Qi; end
        end
        [fitness,idx]=sort(fitness); nidos=nidos(idx,:);
        n_ab = round(pa*n_pop); idx_ab = (n_pop-n_ab+1):n_pop;
        for k = idx_ab, nidos(k,:) = lb + rand(1,dim).*(ub-lb); end
        Qab = zeros(1,numel(idx_ab));
        parfor m = 1:numel(idx_ab), Qab(m) = obj(nidos(idx_ab(m),:)); end %#ok<PFBNS>
        fitness(idx_ab)=Qab; [fitness,idx]=sort(fitness); nidos=nidos(idx,:);
        if fitness(1) < best_fitness, best_fitness=fitness(1); best_nido=nidos(1,:); stall=0;
        else, stall=stall+1; end
        hist_iter(gen)=gen; hist_fval(gen)=best_fitness;
        if mod(gen,10)==0||gen==1, fprintf('  %-4d  %.4f\n', gen, best_fitness); end
        if stall==15, alpha=alpha*2.0; end
        if stall>=30, fprintf('  Parada gen %d\n', gen); break; end
    end
    hist_iter=hist_iter(1:gen); hist_fval=hist_fval(1:gen);
    fprintf('  J_final = %.6f\n\n', best_fitness);
end

function print_gains(t, G, f)
    fprintf('=== %s ===\n', t);
    fprintf('  kpx=%.3f kix=%.4f kdx=%.3f | kpy=%.3f kiy=%.4f kdy=%.3f\n', G.kpx,G.kix,G.kdx,G.kpy,G.kiy,G.kdy);
    fprintf('  kpz=%.3f kiz=%.4f kdz=%.3f | kpphi=%.2f kdphi=%.2f kptheta=%.2f kdtheta=%.2f kppsi=%.2f kdpsi=%.2f\n', ...
        G.kpz,G.kiz,G.kdz,G.kpphi,G.kdphi,G.kptheta,G.kdtheta,G.kppsi,G.kdpsi);
    fprintf('  J = %.6f\n\n', f);
end

function G = vec2gains(v)
    G.kpx=v(1);G.kix=v(2);G.kdx=v(3); G.kpy=v(4);G.kiy=v(5);G.kdy=v(6);
    G.kpz=v(7);G.kiz=v(8);G.kdz=v(9); G.kpphi=v(10);G.kdphi=v(11);
    G.kptheta=v(12);G.kdtheta=v(13); G.kppsi=v(14);G.kdpsi=v(15);
end

function r = rms_(x)   % RMS sin depender de toolboxes
    r = sqrt(mean(x(:).^2));
end

function rd = get_trajectory(t, wp, t_wp, n_seg)
    N=length(t); rd=zeros(N,4);
    for c = 1:4
        for s = 1:n_seg
            T=t_wp(s+1)-t_wp(s);
            if s<n_seg, idx=t>=t_wp(s)&t<t_wp(s+1); else, idx=t>=t_wp(s)&t<=t_wp(s+1); end
            if T<1e-9, rd(idx,c)=wp(s,c); continue; end
            A=[0,0,0,0,0,1; T^5,T^4,T^3,T^2,T,1; 0,0,0,0,1,0;
               5*T^4,4*T^3,3*T^2,2*T,1,0; 0,0,0,2,0,0; 20*T^3,12*T^2,6*T,2,0,0];
            kc=A\[wp(s,c);wp(s+1,c);0;0;0;0]; tau=t(idx)-t_wp(s);
            rd(idx,c)=kc(1)*tau.^5+kc(2)*tau.^4+kc(3)*tau.^3+kc(4)*tau.^2+kc(5)*tau+kc(6);
        end
    end
end

function B = B_segment(r_obs, p1, p2, I)
    mu0=4*pi*1e-7; r_obs=r_obs(:); p1=p1(:); p2=p2(:);
    lv=p2-p1; L=norm(lv); if L<1e-12, B=zeros(3,1); return; end
    lhat=lv/L; r1=r_obs-p1; r2=r_obs-p2; d_vec=r1-(r1'*lhat)*lhat; d=norm(d_vec);
    if d<1e-9, B=zeros(3,1); return; end
    cos1=(r1'*lhat)/norm(r1); cos2=(r2'*lhat)/norm(r2);
    B_mag=(mu0*I)/(4*pi*d)*(cos1-cos2); n_hat=cross(lhat,d_vec/d); B=B_mag*n_hat;
end

function B = B_total(r_obs, wires)
    B=zeros(3,1);
    for s = 1:(size(wires.nodes,1)-1)
        B = B + B_segment(r_obs, wires.nodes(s,:), wires.nodes(s+1,:), wires.I);
    end
end

function [psi_wire, foot] = nearest_segment(r_obs, wires)
    r_obs=r_obs(:); dist_min=inf; psi_wire=0; foot=r_obs;
    for s = 1:(size(wires.nodes,1)-1)
        p1=wires.nodes(s,:)'; p2=wires.nodes(s+1,:)'; lv=p2-p1; L=norm(lv);
        if L<1e-12, continue; end
        t_=max(0,min(1,dot(r_obs-p1,lv)/L^2)); f=p1+t_*lv; d=norm(r_obs-f);
        if d<dist_min, dist_min=d; foot=f; psi_wire=atan2(lv(2),lv(1)); end
    end
end

function precomp = make_precomp(t, rd, mag, wires, Ts)
    N=length(t); mag_step=max(1,round(mag.Ts/Ts));
    tau_np=1.0; a_np=mag.Ts/(tau_np+mag.Ts); nx_f=[]; ny_f=[];
    precomp.psi_wire=zeros(N,1); precomp.B_ref_vec=zeros(N,3);
    precomp.dBdperp=zeros(N,3); precomp.dBdz=zeros(N,3); precomp.dthe_dperp=zeros(N,1);
    last_psi=0; last_Bref=zeros(3,1); last_dBp=zeros(3,1); last_dBz=zeros(3,1); last_tp=0;
    for k = 1:N
        if mod(k-1,mag_step)==0
            r_ref=[rd(k,1);rd(k,2);rd(k,3)];
            [psi_w,~]=nearest_segment(r_ref,wires);
            nx_raw=-sin(psi_w); ny_raw=cos(psi_w);
            if isempty(nx_f), nx_f=nx_raw; ny_f=ny_raw;
            else, nx_f=nx_f+a_np*(nx_raw-nx_f); ny_f=ny_f+a_np*(ny_raw-ny_f); end
            nn=hypot(nx_f,ny_f); psi_w=atan2(-nx_f/nn, ny_f/nn);
            n_perp=[-sin(psi_w);cos(psi_w);0]; n_z=[0;0;1];
            Bref=B_total(r_ref,wires);
            Bp=B_total(r_ref+mag.dr*n_perp,wires); Bm=B_total(r_ref-mag.dr*n_perp,wires);
            Bzp=B_total(r_ref+mag.dr*n_z,wires);   Bzm=B_total(r_ref-mag.dr*n_z,wires);
            wrp=@(a) mod(a+pi,2*pi)-pi;
            thp=atan2(Bp(2),Bp(1)); thm=atan2(Bm(2),Bm(1));
            last_tp=wrp(thp-thm)/(2*mag.dr);            % dtheta/dperp [rad/m]
            last_psi=psi_w; last_Bref=Bref;
            last_dBp=(Bp-Bm)/(2*mag.dr); last_dBz=(Bzp-Bzm)/(2*mag.dr);
        end
        precomp.psi_wire(k)=last_psi; precomp.B_ref_vec(k,:)=last_Bref';
        precomp.dBdperp(k,:)=last_dBp'; precomp.dBdz(k,:)=last_dBz'; precomp.dthe_dperp(k)=last_tp;
    end
end

function [om,Im,Vm,Pm,Pe] = actuator_model(F,t2,t3,t4,xd_,yd_,zd_,phi,th,psi,dphi,dth,dps,p)
    w2=max(p.Mact_inv*[F;t2;t3;t4],0); om=sqrt(w2); Im=p.kdrot*w2/p.km; Vm=p.R*Im+p.km*om;
    Rez=[cos(psi)*sin(th)+sin(psi)*sin(phi)*cos(th); sin(psi)*sin(th)-cos(psi)*sin(phi)*cos(th); cos(phi)*cos(th)];
    Pm=F*(Rez'*[xd_;yd_;zd_])+t2*dphi+t3*dth+t4*dps; Pe=p.R*sum(Im.^2);
end

function dx = sys_deriv(x,F,t2,t3,t4,psid,mp,Ix,Iy,Iz,p,pert,vf)
    xdot=x(2);ydot=x(4);zdot=x(6); phi=x(7);phidot=x(8); th=x(9);thdot=x(10); psidot=x(12);
    xdd=p.g*(th*cos(psid)+phi*sin(psid)); ydd=p.g*(th*sin(psid)-phi*cos(psid)); zdd=F/mp-p.g;
    phidd=t2/Ix; thdd=t3/Iy; psidd=t4/Iz;
    if pert
        ve=[xdot;ydot;zdot]+0.05*[phidot;thdot;0]-vf;
        Fa1=-p.k1e*sqrt(ve(1)^2+ve(2)^2)*ve(1)-p.k2e*abs(ve(1))*ve(1);
        Fa2=-p.k1e*sqrt(ve(1)^2+ve(2)^2)*ve(2)-p.k2e*abs(ve(2))*ve(2);
        Fa3=-p.k3e*(ve(1)^2+ve(2)^2)-p.k4e*abs(ve(3))*ve(3);
        Fae=[Fa1;Fa2;Fa3]-0.28*sqrt(abs(F))*vf; Mae=-0.05*cross([0;0;1],Fae);
        xdd=xdd+Fae(1)/mp; ydd=ydd+Fae(2)/mp; zdd=zdd+Fae(3)/mp;
        phidd=phidd+Mae(1)/Ix; thdd=thdd+Mae(2)/Iy; psidd=psidd+Mae(3)/Iz;
    end
    dx=[xdot;xdd;ydot;ydd;zdot;zdd;phidot;phidd;thdot;thdd;psidot;psidd];
end

%% ── mag_solve_perp: estimador de la desviacion PERP, grad o ang ─────
function [delta, w] = mag_solve_perp(B_meas, k, precomp, mag, method)
    B_meas=B_meas(:); B_ref=precomp.B_ref_vec(k,:)';
    if strcmp(method,'grad')
        Jp=precomp.dBdperp(k,:)'; np2=Jp'*Jp; lam=mag.sigB^2/mag.sigma_s^2;
        e_B=B_ref-B_meas; delta=-(Jp'*e_B)/(np2+lam); w=np2/(np2+lam);
    else % 'ang': solo direccion del campo
        wrp=@(a) mod(a+pi,2*pi)-pi;
        th_meas=atan2(B_meas(2),B_meas(1)); th_ref=atan2(B_ref(2),B_ref(1));
        r_th=wrp(th_meas-th_ref); ga=precomp.dthe_dperp(k);
        sig_th=mag.sigB/max(norm(B_ref),1e-12); lam_a=sig_th^2/mag.sigma_s^2;
        delta=(ga*r_th)/(ga^2+lam_a); w=ga^2/(ga^2+lam_a);
    end
    delta=max(min(delta,5),-5);
end

function [delta, w] = mag_estimate(r_obs, k, precomp, mag, wires, rs, method)
    B_meas=B_total(r_obs,wires)+mag.sigB*randn(rs,3,1);
    [delta, w]=mag_solve_perp(B_meas, k, precomp, mag, method);
end

%% ── run_sim_fast: PLANTA LIMPIA (perp sin deriva) para optimizar ────
function [J,div] = run_sim_fast(G,p,t,rd,Ts,N,x0,err_max,precomp,~,~,~,~)
    kpx=G.kpx;kix=G.kix;kdx=G.kdx; kpy=G.kpy;kiy=G.kiy;kdy=G.kdy; kpz=G.kpz;kiz=G.kiz;kdz=G.kdz;
    kph=G.kpphi;kdph=G.kdphi; kth=G.kptheta;kdth=G.kdtheta; kps=G.kppsi;kdps=G.kdpsi;
    g_=p.g; mg_=p.m*p.g; m_=p.m; Ix_=p.Ixx;Iy_=p.Iyy;Iz_=p.Izz; Mi=p.Mact_inv;kdr=p.kdrot;km_=p.km;R_=p.R;
    pi4=deg2rad(40); Ts2=Ts*.5; Ts6=Ts/6; ilim=100;
    x=x0; int_ex=0;int_ey=0;int_ez=0; Jacc=0; div=false; J=0;
    for k=1:N
        xd=rd(k,1);yd=rd(k,2);zd=rd(k,3);psd=rd(k,4);
        xpos=x(1);ypos=x(3);zpos=x(5); phi=x(7);th=x(9);psi=x(11);
        ex=xd-xpos; ey=yd-ypos; ez=zd-zpos;
        if abs(ex)>err_max||abs(ey)>err_max||abs(ez)>err_max||abs(phi)>1.5||abs(th)>1.5||~isfinite(x(1)), div=true; return; end
        int_ex=max(min(int_ex+ex*Ts,ilim),-ilim); int_ey=max(min(int_ey+ey*Ts,ilim),-ilim); int_ez=max(min(int_ez+ez*Ts,ilim),-ilim);
        vx=x(2);vy=x(4);vz=x(6);dphi=x(8);dth=x(10);dps=x(12);
        F=max(mg_+kpz*ez+kiz*int_ez-kdz*vz,0.05);
        thd=min(max(kpx*ex+kix*int_ex-kdx*vx,-pi4),pi4);
        phd=min(max(-(kpy*ey+kiy*int_ey-kdy*vy),-pi4),pi4);
        ephi=phd-phi;eth=thd-th;eps_=psd-psi;
        t2=kph*ephi-kdph*dphi; t3=kth*eth-kdth*dth; t4=kps*eps_-kdps*dps;
        w2=max(Mi*[F;t2;t3;t4],0); Im2=kdr*w2/km_;
        Rez=[cos(psi)*sin(th)+sin(psi)*sin(phi)*cos(th);sin(psi)*sin(th)-cos(psi)*sin(phi)*cos(th);cos(phi)*cos(th)];
        Pm=F*(Rez'*[vx;vy;vz])+t2*dphi+t3*dth+t4*dps; Pe=R_*sum(Im2.^2);
        Jacc=Jacc+Ts*(t(k)*(ex^2+ey^2+ez^2+ephi^2+eth^2+eps_^2)+p.wPm*Pm^2+p.wPe*Pe+p.wRate*(dphi^2+dth^2));
        cp=cos(psd);sp=sin(psd); Fmg=F/m_-g_; t2i=t2/Ix_;t3i=t3/Iy_;t4i=t4/Iz_;
        f1=[x(2);g_*(x(9)*cp+x(7)*sp);x(4);g_*(x(9)*sp-x(7)*cp);x(6);Fmg;x(8);t2i;x(10);t3i;x(12);t4i];
        x2=x+Ts2*f1; f2=[x2(2);g_*(x2(9)*cp+x2(7)*sp);x2(4);g_*(x2(9)*sp-x2(7)*cp);x2(6);Fmg;x2(8);t2i;x2(10);t3i;x2(12);t4i];
        x3=x+Ts2*f2; f3=[x3(2);g_*(x3(9)*cp+x3(7)*sp);x3(4);g_*(x3(9)*sp-x3(7)*cp);x3(6);Fmg;x3(8);t2i;x3(10);t3i;x3(12);t4i];
        x4=x+Ts*f3;  f4=[x4(2);g_*(x4(9)*cp+x4(7)*sp);x4(4);g_*(x4(9)*sp-x4(7)*cp);x4(6);Fmg;x4(8);t2i;x4(10);t3i;x4(12);t4i];
        x=x+Ts6*(f1+2*f2+2*f3+f4);
    end
    J=Jacc;
end

%% ── run_sim: 1D perp. method = 'none' | 'grad' | 'ang' ──────────────
function S = run_sim(G,p,t,rd,Ts,N,x0,precomp,mag,wires,nav,method)
    mp=p.m;Ixxp=p.Ixx;Iyyp=p.Iyy;Izzp=p.Izz;vf=[0;0;0];   % sin perturbaciones
    int_lim=100; x_st=x0; mag_step=max(1,round(mag.Ts/Ts));
    use_mag = ~strcmp(method,'none');
    if use_mag, rs=RandStream('mt19937ar','Seed',7); end
    X=zeros(N,12);U=zeros(N,4);Ang=zeros(N,2);Jo=zeros(N,1);Pmo=zeros(N,1);Peo=zeros(N,1);
    dlt=zeros(N,1);etap=zeros(N,1);wp_=zeros(N,1);dpt=zeros(N,1);
    int_ex=0;int_ey=0;int_ez=0;J_acc=0; delta=0;w_perp=0;mag_ctr=0;eta_perp=0; S.diverged=false;
    for k=1:N
        xd=rd(k,1);yd=rd(k,2);zd=rd(k,3);psid=rd(k,4);
        xpos=x_st(1);xdot=x_st(2);ypos=x_st(3);ydot=x_st(4);zpos=x_st(5);zdot=x_st(6);
        phi=x_st(7);phidot=x_st(8);theta=x_st(9);thetadot=x_st(10);psi=x_st(11);psidot=x_st(12);
        psi_wire=precomp.psi_wire(k); spw=sin(psi_wire);cpw=cos(psi_wire);
        v_drift=0.5*p.g*nav.b_g*t(k)^2;
        dp_true=(xpos-xd)*(-spw)+(ypos-yd)*cpw;
        if use_mag
            if mag_ctr==0, [delta,w_perp]=mag_estimate([xpos;ypos;zpos],k,precomp,mag,wires,rs,method); end
            mag_ctr=mod(mag_ctr+1,mag_step);
            d_hat=dp_true+eta_perp;
            eta_perp=eta_perp+Ts*(v_drift-nav.k_obs*w_perp*(d_hat-delta));
        else
            eta_perp=eta_perp+Ts*v_drift; w_perp=0;
        end
        ex=(xd-xpos)+eta_perp*spw; ey=(yd-ypos)-eta_perp*cpw; ez=zd-zpos;
        exd=-xdot;eyd=-ydot;ezd=-zdot;
        dlt(k)=delta; etap(k)=eta_perp; wp_(k)=w_perp; dpt(k)=dp_true;
        int_ex=max(min(int_ex+ex*Ts,int_lim),-int_lim); int_ey=max(min(int_ey+ey*Ts,int_lim),-int_lim); int_ez=max(min(int_ez+ez*Ts,int_lim),-int_lim);
        F=max(p.m*p.g+G.kpz*ez+G.kiz*int_ez+G.kdz*ezd,0.05);
        theta_d=min(max(G.kpx*ex+G.kix*int_ex+G.kdx*exd,-deg2rad(40)),deg2rad(40));
        phi_d=min(max(-(G.kpy*ey+G.kiy*int_ey+G.kdy*eyd),-deg2rad(40)),deg2rad(40));
        e_phi=phi_d-phi;e_theta=theta_d-theta;e_psi=psid-psi;
        tau2=G.kpphi*e_phi-G.kdphi*phidot; tau3=G.kptheta*e_theta-G.kdtheta*thetadot; tau4=G.kppsi*e_psi-G.kdpsi*psidot;
        U(k,:)=[F,tau2,tau3,tau4]; Ang(k,:)=[phi_d,theta_d];
        [~,~,~,Pm,Pe]=actuator_model(F,tau2,tau3,tau4,xdot,ydot,zdot,phi,theta,psi,phidot,thetadot,psidot,p);
        Pmo(k)=Pm; Peo(k)=Pe;
        f1=sys_deriv(x_st,F,tau2,tau3,tau4,psid,mp,Ixxp,Iyyp,Izzp,p,false,vf);
        f2=sys_deriv(x_st+Ts/2*f1,F,tau2,tau3,tau4,psid,mp,Ixxp,Iyyp,Izzp,p,false,vf);
        f3=sys_deriv(x_st+Ts/2*f2,F,tau2,tau3,tau4,psid,mp,Ixxp,Iyyp,Izzp,p,false,vf);
        f4=sys_deriv(x_st+Ts*f3,F,tau2,tau3,tau4,psid,mp,Ixxp,Iyyp,Izzp,p,false,vf);
        x_st=x_st+(Ts/6)*(f1+2*f2+2*f3+f4); X(k,:)=x_st';
        exr=xd-xpos;eyr=yd-ypos;ezr=zd-zpos; r_t=[exr;eyr;ezr;e_phi;e_theta;e_psi];
        J_acc=J_acc+Ts*(t(k)*(r_t'*r_t)+p.wPm*Pm^2+p.wPe*Pe+p.wRate*(phidot^2+thetadot^2)); Jo(k)=J_acc;
    end
    S.X=X;S.U=U;S.Ang=Ang;S.J=Jo;S.Pm=Pmo;S.Pe=Peo;S.delta=dlt;S.eta_perp=etap;S.w_perp=wp_;S.dp_true=dpt;
end

function J = obj_quad(v,p,rd,t,Ts,N,x0,err_max,precomp,mag,wires,nav,~)
    G=vec2gains(v); [J,div]=run_sim_fast(G,p,t,rd,Ts,N,x0,err_max,precomp,mag,wires,nav,'none');
    if div||~isfinite(J), J=1e8; end
end

function plot_convergencia(hi, hf, fval, G, metodo, fp)
    if numel(hf)<2, return; end
    [vm,im]=min(hf);
    fig=figure('Position',[100 80 1000 460],'Color','w');
    semilogy(hi,hf,'b-o','MarkerSize',3,'LineWidth',1.8); hold on; semilogy(im,vm,'r*','MarkerSize',14,'LineWidth',2.5);
    xlabel('Iteracion'); ylabel('J (log)'); title(['Convergencia — ' metodo]); grid on;
    str=sprintf('kpx=%.2f kdx=%.2f\nkpy=%.2f kdy=%.2f\nkpz=%.2f kdz=%.2f\nkpphi=%.1f kdphi=%.1f\nkptheta=%.1f kdtheta=%.1f\nJ=%.3f',...
        G.kpx,G.kdx,G.kpy,G.kdy,G.kpz,G.kdz,G.kpphi,G.kdphi,G.kptheta,G.kdtheta,fval);
    annotation('textbox',[0.68 0.5 0.2 0.3],'String',str,'FitBoxToText','on','BackgroundColor',[1 1 .88],'FontName','Courier New','FontSize',8,'Interpreter','none');
    drawnow; try exportgraphics(fig,fp,'Resolution',300); catch, print(fig,fp,'-dpng','-r300'); end; close(fig);
end

function plot_results(S, t, rd, Tfin, titulo, fp)
    x=S.X(:,1);y=S.X(:,3);z=S.X(:,5); phi=S.X(:,7);theta=S.X(:,9);psi=S.X(:,11);
    phi_d=S.Ang(:,1);th_d=S.Ang(:,2); ex=rd(:,1)-x;ey=rd(:,2)-y;ez=rd(:,3)-z;
    ephi=phi_d-phi;eth=th_d-theta;epsi=rd(:,4)-psi; F=S.U(:,1);tau2=S.U(:,2);tau3=S.U(:,3);tau4=S.U(:,4); xl=[0 Tfin];
    fig=figure('Position',[20 20 1600 1200],'Color','w','Name',titulo);
    subplot(5,3,1);plot(t,rd(:,1),'r--',t,x,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('x[m]');title('Posicion X');   % leyenda quitada
    subplot(5,3,2);plot(t,rd(:,2),'r--',t,y,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('y[m]');title('Posicion Y');
    subplot(5,3,3);plot(t,rd(:,3),'r--',t,z,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('z[m]');title('Posicion Z');
    subplot(5,3,4);plot(t,ex,'b','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('e_x[m]');title('Error e_x');
    subplot(5,3,5);plot(t,ey,'r','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('e_y[m]');title('Error e_y');
    subplot(5,3,6);plot(t,ez,'g','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('e_z[m]');title('Error e_z');
    subplot(5,3,7);plot(t,phi_d*180/pi,'r--',t,phi*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('phi[deg]');title('Alabeo');   % leyenda quitada
    subplot(5,3,8);plot(t,th_d*180/pi,'r--',t,theta*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('theta[deg]');title('Cabeceo');
    subplot(5,3,9);plot(t,rd(:,4)*180/pi,'r--',t,psi*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);ylabel('psi[deg]');title('Guinada');
    subplot(5,3,10);plot(t,ephi*180/pi,'b','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('e_phi[deg]');title('Error alabeo');
    subplot(5,3,11);plot(t,eth*180/pi,'r','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('e_theta[deg]');title('Error cabeceo');
    subplot(5,3,12);plot(t,epsi*180/pi,'g','LineWidth',1.2);grid on;xlim(xl);yline(0,'k:');ylabel('e_psi[deg]');title('Error guinada');
    subplot(5,3,13);yyaxis left;plot(t,F,'b','LineWidth',1.2);ylabel('F[N]');yyaxis right;plot(t,S.Pm,'r--',t,S.Pe,'g--');ylabel('P[W]');grid on;xlim(xl);xlabel('t[s]');title('Empuje y Potencias');legend('F','Pm','Pe');
    subplot(5,3,14);plot(t,tau2,'b',t,tau3,'r',t,tau4,'g','LineWidth',1.2);grid on;xlim(xl);ylabel('[Nm]');xlabel('t[s]');legend('tau2','tau3','tau4');title('Torques');
    subplot(5,3,15);plot(t,S.J,'m','LineWidth',1.5);grid on;xlim(xl);ylabel('J');xlabel('t[s]');title(sprintf('J_{fin}=%.3f',S.J(end)));
    set(findall(fig,'-property','FontSize'),'FontSize',8); sgtitle(titulo,'FontSize',11,'FontWeight','bold'); drawnow;
    try exportgraphics(fig,fp,'Resolution',300); catch, print(fig,fp,'-dpng','-r300'); end; close(fig);
end

%% ── plot_compara_1D: la grafica clave de la comparacion ─────────────
function plot_compara_1D(Sc, Sg, Sa, t, fp)
    fig=figure('Position',[60 60 1200 800],'Color','w','Name','Comparativa 1D grad vs ang');
    subplot(2,1,1);
    plot(t, Sc.dp_true,'k','LineWidth',1.4); hold on;
    plot(t, Sg.dp_true,'b','LineWidth',1.6);
    plot(t, Sa.dp_true,'r','LineWidth',1.6); grid on; yline(0,'k:');
    ylabel('desviacion perp REAL [m]'); xlabel('t [s]');
    title(sprintf('Desviacion perpendicular real  |  RMS: control=%.3f  grad=%.3f  ang=%.3f m', ...
        rms_(Sc.dp_true), rms_(Sg.dp_true), rms_(Sa.dp_true)));
    legend('control (sin correccion)','GRAD','ANG','Location','best');
    subplot(2,1,2);
    plot(t, Sc.J,'k','LineWidth',1.4); hold on;
    plot(t, Sg.J,'b','LineWidth',1.6);
    plot(t, Sa.J,'r','LineWidth',1.6); grid on;
    ylabel('Indice J acumulado'); xlabel('t [s]');
    title(sprintf('Coste J  |  fin: control=%.1f  grad=%.1f  ang=%.1f', Sc.J(end), Sg.J(end), Sa.J(end)));
    legend('control','GRAD','ANG','Location','best');
    sgtitle('Comparativa correccion PERP: gradiente vs angulo (L, 1D, sin pert)','FontSize',12,'FontWeight','bold');
    drawnow; try exportgraphics(fig,fp,'Resolution',300); catch, print(fig,fp,'-dpng','-r300'); end; close(fig);
end
