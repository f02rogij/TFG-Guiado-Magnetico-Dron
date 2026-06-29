%% ================================================================
%  campos3.m   (nueva version de campos2.m)
%  Comparacion de TRES metodos para calcular el campo magnetico B, en
%  TIEMPO y en CAMPO, sobre DOS distribuciones de corriente:
%     - escalera en esquina   (railes antiparalelos en L + travesanos)
%     - L simple              (un unico hilo en L)
%
%  Metodos:
%   (1) ANALITICO : formula cerrada del segmento recto finito (la de los
%                   dron_*.m). EXACTO para hilos rectos -> es la REFERENCIA.
%   (2) NUMERICO  : integracion directa de Biot-Savart discretizando los
%                   hilos en trozos de longitud ds (la de campos2.m).
%   (3) TABLA     : se precalcula B en una malla 3D (con el analitico) y en
%                   cada consulta se INTERPOLA (CUBICA: gradiente suave y
%                   menor error que trilineal). O(1) por llamada,
%                   independiente del nº de conductores -> ideal para las
%                   miles de llamadas del dron.
%
%  Salida: una figura por distribucion en la carpeta  tfg\campo\
%% ================================================================
clear; close all;
[sd, ~] = fileparts(mfilename('fullpath'));
if isempty(sd), sd = pwd; end
outdir = fullfile(sd, 'campo');                 % <<< carpeta campo
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('Salida en: %s\n', outdir);

%% ── Parametros de los metodos ───────────────────────────────────────
ds = 0.10;     % paso de discretizacion del NUMERICO [m]
h  = 0.50;     % paso de la malla de la TABLA [m]

%% ── Distribucion 1: ESCALERA EN ESQUINA ────────────────────────────
Lc=30; w=4; P=8; r1x=3:P:21; r2y=9:P:Lc;
cond={};
cond{end+1}=struct('nodes',[-2 -w 0; Lc+w -w 0; Lc+w Lc+2 0],'sig',+1);   % riel exterior
cond{end+1}=struct('nodes',[-2 +w 0; Lc-w +w 0; Lc-w Lc+2 0],'sig',-1);   % riel interior
for k=1:numel(r1x), cond{end+1}=struct('nodes',[r1x(k) -w 0; r1x(k) +w 0],'sig',(-1)^k); end %#ok<SAGROW>
for k=1:numel(r2y), cond{end+1}=struct('nodes',[Lc-w r2y(k) 0; Lc+w r2y(k) 0],'sig',(-1)^(k+1)); end %#ok<SAGROW>
w_esc.cond=cond; w_esc.I=200;

%% ── Distribucion 2: L SIMPLE (un hilo) ─────────────────────────────
Lsz=15;
w_L.cond = { struct('nodes',[0 0 0; Lsz 0 0; Lsz Lsz 0],'sig',+1) };
w_L.I    = 200;

%% ── Lista de casos (geometria + dominio de la tabla + linea de test) ─
casos = {};
% La linea de corte se mantiene DENTRO del dominio de la tabla (si no, la
% interpolacion extrapola en los bordes y dispara un error artificial).
casos{end+1} = struct('nombre','escalera_esquina', 'wires',w_esc, ...
    'dom', struct('x',[-2 33],'y',[-6 34],'z',[3 7]), ...
    'linea', struct('x',10,  'z',5, 'y',linspace(-6,14,600)'));
casos{end+1} = struct('nombre','L_simple', 'wires',w_L, ...
    'dom', struct('x',[-5 20],'y',[-5 20],'z',[3 7]), ...
    'linea', struct('x',7.5, 'z',5, 'y',linspace(-5,10,600)'));

%% ── Ejecutar la comparacion en cada distribucion ────────────────────
for ic = 1:numel(casos)
    comparar_caso(casos{ic}, ds, h, outdir);
end
fprintf('\nCompletado.\n');


%% ================================================================
%%  FUNCIONES LOCALES
%% ================================================================

%% ── comparar_caso: hace toda la comparacion de una distribucion ────
function comparar_caso(caso, ds, h, outdir)
    wires=caso.wires; dom=caso.dom; nombre=caso.nombre; n_cond=numel(wires.cond);
    disp_n=strrep(nombre,'_',' ');   % nombre para titulos (sin '_' -> sin subindice TeX)
    fprintf('\n===== DISTRIBUCION: %s (%d conductores, I=%d A) =====\n', nombre, n_cond, wires.I);

    % (3) construir tabla (coste de setup)
    tic; [Fx,Fy,Fz,n_nodos]=build_tabla(wires,dom,h); t_build=toc;
    fprintf('  tabla: %d nodos en %.3f s\n', n_nodos, t_build);

    % puntos de consulta para el benchmark
    Nq=1000; rng(1);
    Q=[dom.x(1)+diff(dom.x)*rand(Nq,1), dom.y(1)+diff(dom.y)*rand(Nq,1), dom.z(1)+diff(dom.z)*rand(Nq,1)];
    Ba=zeros(Nq,3); Bn=zeros(Nq,3); Bt=zeros(Nq,3);
    tic; for i=1:Nq, Ba(i,:)=B_analitico(Q(i,:),wires)'; end; t_ana=toc;
    tic; for i=1:Nq, Bn(i,:)=B_numerico (Q(i,:),wires,ds)'; end; t_num=toc;
    tic; for i=1:Nq, Bt(i,:)=[Fx(Q(i,1),Q(i,2),Q(i,3)),Fy(Q(i,1),Q(i,2),Q(i,3)),Fz(Q(i,1),Q(i,2),Q(i,3))]; end; t_tab=toc;

    us=1e6; uT=1e6;
    fprintf('  t/llamada [us]:  analitico=%.2f  numerico=%.2f  tabla=%.2f  (setup %.3f s)\n', ...
        us*t_ana/Nq, us*t_num/Nq, us*t_tab/Nq, t_build);

    % errores vs analitico (referencia exacta)
    en=sqrt(sum((Bn-Ba).^2,2))*uT; et=sqrt(sum((Bt-Ba).^2,2))*uT; Bmag=sqrt(sum(Ba.^2,2))*uT;
    fprintf('  error RMS vs analitico:  numerico=%.4f uT   tabla=%.4f uT  (|B| medio=%.2f uT)\n', ...
        rms_(en), rms_(et), mean(Bmag));

    % cruce de coste tiempo total vs numero de llamadas K
    K=round(logspace(1,7,60)); tot_ana=K*(t_ana/Nq); tot_num=K*(t_num/Nq); tot_tab=t_build+K*(t_tab/Nq);
    dif=(t_ana-t_tab)/Nq;                       % ahorro por llamada de la tabla vs analitico
    if dif>0, kcross=t_build/dif; else, kcross=Inf; end   % Inf = la tabla no supera al analitico

    % linea de prueba (corte perpendicular)
    yl=caso.linea.y; xl=caso.linea.x; zl=caso.linea.z;
    La=zeros(numel(yl),3); Ln=zeros(numel(yl),3); Lt=zeros(numel(yl),3);
    for i=1:numel(yl)
        La(i,:)=B_analitico([xl yl(i) zl],wires)';
        Ln(i,:)=B_numerico ([xl yl(i) zl],wires,ds)';
        Lt(i,:)=[Fx(xl,yl(i),zl), Fy(xl,yl(i),zl), Fz(xl,yl(i),zl)];
    end
    mg=@(M) sqrt(sum(M.^2,2))*uT;

    % ── figura ──
    fig=figure('Position',[40 40 1500 900],'Color','w','Name',['campos3 — ' nombre]);
    subplot(2,3,1);
    bar([us*t_ana/Nq, us*t_num/Nq, us*t_tab/Nq]); grid on; set(gca,'XTickLabel',{'Analitico','Numerico','Tabla'});
    ylabel('t por llamada [\mus]'); title('Coste por consulta'); set(gca,'YScale','log');

    subplot(2,3,2);
    loglog(K,tot_ana,'b','LineWidth',1.8); hold on; loglog(K,tot_num,'r','LineWidth',1.8); loglog(K,tot_tab,'g','LineWidth',1.8);
    allt=[tot_ana tot_num tot_tab]; xlim([K(1) K(end)]); ylim([min(allt) max(allt)]);   % ejes al rango de datos
    if isfinite(kcross) && kcross<=max(K)
        xline(kcross,'k:','LineWidth',1.2);
        tit_ct=sprintf('Coste total (tabla<analitico si K>%.0f)',kcross);
    else
        tit_ct='Coste total (tabla no supera al analitico)';
    end
    grid on; xlabel('numero de llamadas K'); ylabel('tiempo total [s]'); title(tit_ct);
    legend('Analitico','Numerico','Tabla','Location','northwest');

    subplot(2,3,3);
    histogram(en,'FaceColor','r','FaceAlpha',0.5); hold on; histogram(et,'FaceColor','g','FaceAlpha',0.5);
    grid on; xlabel('error vs analitico [\muT]'); ylabel('nº consultas'); title('Distribucion del error'); legend('Numerico','Tabla');

    subplot(2,3,4);
    plot(yl,mg(La),'b','LineWidth',2); hold on; plot(yl,mg(Ln),'r--','LineWidth',1.2); plot(yl,mg(Lt),'g-.','LineWidth',1.2);
    grid on; xlim([yl(1) yl(end)]); xlabel(sprintf('y [m]  (corte x=%g, z=%g)',xl,zl)); ylabel('|B| [\muT]'); title('Campo en la linea de prueba');
    legend('Analitico','Numerico','Tabla','Location','best');

    subplot(2,3,5);
    plot(yl, mg(Ln-La),'r','LineWidth',1.4); hold on; plot(yl, mg(Lt-La),'g','LineWidth',1.4);
    grid on; xlim([yl(1) yl(end)]); xlabel('y [m]'); ylabel('|error| [\muT]'); title('Error en la linea (vs analitico)'); legend('Numerico','Tabla','Location','best');

    if isfinite(kcross), s_cross=sprintf('Tabla gana a analitico si K>%.0f',kcross);
    else, s_cross='Tabla NO supera al analitico (geom. simple)'; end
    subplot(2,3,6); axis off;
    txt=sprintf(['DISTRIBUCION: %s\n\nConductores: %d   I=%d A\nnumerico ds=%.2f m\ntabla h=%.2f m (%d nodos)\n\n' ...
        't/llamada [us]:\n  Analitico = %.2f\n  Numerico  = %.2f\n  Tabla     = %.2f\n\nSetup tabla = %.3f s\n%s\n\n' ...
        'Error RMS vs analitico:\n  Numerico = %.4f uT\n  Tabla    = %.4f uT'], ...
        disp_n, n_cond, wires.I, ds, h, n_nodos, us*t_ana/Nq, us*t_num/Nq, us*t_tab/Nq, t_build, s_cross, rms_(en), rms_(et));
    text(0.02,0.98,txt,'VerticalAlignment','top','FontName','Courier New','FontSize',9,'Interpreter','none');

    sgtitle(sprintf('campos3 — %s : Analitico vs Numerico vs Tabla',disp_n),'FontSize',13,'FontWeight','bold');
    drawnow;
    fp=fullfile(outdir, sprintf('comparativa_metodos_%s.png', nombre));
    try exportgraphics(fig, fp,'Resolution',300); catch, print(fig, fp,'-dpng','-r300'); end
    close(fig);
    fprintf('  figura: %s\n', fp);
end

%% ── (1) ANALITICO: segmento recto finito (formula cerrada) ─────────
function B = B_analitico(r_obs, wires)
    B=zeros(3,1);
    for c=1:numel(wires.cond)
        nd=wires.cond{c}.nodes; Ic=wires.cond{c}.sig*wires.I;
        for s=1:size(nd,1)-1
            B=B+seg_analitico(r_obs, nd(s,:), nd(s+1,:), Ic);
        end
    end
end
function B = seg_analitico(r_obs, p1, p2, I)
    mu0=4*pi*1e-7; r_obs=r_obs(:); p1=p1(:); p2=p2(:);
    lv=p2-p1; L=norm(lv); if L<1e-12, B=zeros(3,1); return; end
    lhat=lv/L; r1=r_obs-p1; r2=r_obs-p2; d_vec=r1-(r1'*lhat)*lhat; d=norm(d_vec);
    if d<1e-9, B=zeros(3,1); return; end
    cos1=(r1'*lhat)/norm(r1); cos2=(r2'*lhat)/norm(r2);
    B=(mu0*I)/(4*pi*d)*(cos1-cos2)*cross(lhat,d_vec/d);
end

%% ── (2) NUMERICO: integracion directa de Biot-Savart ──────────────
function B = B_numerico(r_obs, wires, ds)
    mu0=4*pi*1e-7; r_obs=r_obs(:); B=zeros(3,1);
    for c=1:numel(wires.cond)
        nd=wires.cond{c}.nodes; Ic=wires.cond{c}.sig*wires.I;
        for s=1:size(nd,1)-1
            p1=nd(s,:)'; p2=nd(s+1,:)'; seg=p2-p1; L=norm(seg);
            n=max(1,ceil(L/ds)); dl=seg/n;
            for i=1:n
                mid=p1+(i-0.5)*dl; rp=r_obs-mid; rn=norm(rp);
                if rn<1e-9, continue; end
                B=B+(mu0*Ic/(4*pi))*cross(dl,rp)/rn^3;
            end
        end
    end
end

%% ── (3) TABLA: precalculo en malla + interpolantes CUBICOS ─────────
function [Fx,Fy,Fz,n_nodos] = build_tabla(wires, dom, h)
    gx=dom.x(1):h:dom.x(2); gy=dom.y(1):h:dom.y(2); gz=dom.z(1):h:dom.z(2);
    [GX,GY,GZ]=ndgrid(gx,gy,gz);
    Bx=zeros(size(GX)); By=Bx; Bz=Bx;
    for idx=1:numel(GX)
        B=B_analitico([GX(idx);GY(idx);GZ(idx)], wires);   % generador exacto
        Bx(idx)=B(1); By(idx)=B(2); Bz(idx)=B(3);
    end
    Fx=griddedInterpolant(GX,GY,GZ,Bx,'cubic');
    Fy=griddedInterpolant(GX,GY,GZ,By,'cubic');
    Fz=griddedInterpolant(GX,GY,GZ,Bz,'cubic');
    n_nodos=numel(GX);
end

function r = rms_(x), r=sqrt(mean(x(:).^2)); end
