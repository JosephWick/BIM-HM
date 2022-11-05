% BP1_visco_hm2.m

function varargout = bp1vhm (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% XXX ------- Public Functions --------- XXX

% set up h-matrices/kernels
function r = build()
  % this function defines the physical ss.probLem as well as the key value files
  % that hmmvp will use

  addpaths();

  % these are QDBIM style defns
  % shear wave speed (m/s)
  Vs = 3464;
  % density (kg/m^3)
  rho = 2670;
  % shear modulus (MPa)
  %G = rho*Vs^2/1e6;
  G = 30e3; ss.G = G;

  %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                                                      %
  %                        M E S H                       %
  %                                                      %
  % % % % % % % % % % % % % % % % % % % % % % % % % % %  %
  %
  % Naming convention:
  %   - '_c' suffixed terms refer to the center of patches
  %   - 'hat' suffixed terms refer to a 1D vector that represents a direction
  %         of a 2D or 3D mesh
  %

  disp('begin mesh...')

  % Fault Meshes
  ss.probL = 200e3;
  ss.probW = 200e3;


  ss.lambdaZ = 40e3; % fault depth extent
  ss.M = 400; %number of fault cells
  ss.dz = ss.lambdaZ/ss.M; dz = ss.dz;

  ss.transition = 40e3; Transition = ss.transition;
  ss.Ny = 51;
  ss.Nz = 51;

  % FAULT
  % fault patch edges (top left)
  yf = 0;
  faultX = zeros(1,ss.M);
  faultY = zeros(1,ss.M);
  faultZ = linspace(0, ss.lambdaZ-ss.dz, ss.M); faultZ = faultZ';

  % fault patch centers
  faultX_c = faultX;
  faultY_c = faultY;
  faultZ_c = faultZ'+(ss.dz/2);

  % width of fault patches
  Lf = zeros(ss.M,1);
  Wf = ones(ss.M,1)*ss.dz;

  % SHEAR
  % *hat terms are 1xN
  eps = 1e-12;
  nc = (-ss.Ny/2:ss.Ny/2);
  shearZhat = ss.transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*ss.transition;
  shearYhat = tan(nc*pi/(2.5*max(nc)))*32e3;
  shearXhat = zeros(size(shearZhat));

  % shear patch centers
  shearX_c = zeros(ss.Ny,ss.Nz);
  ss.shearY_chat = zeros(ss.Ny,1);
  ss.shearZ_chat = zeros(ss.Nz,1);

  % sizes of shear patches
  L = zeros(1,ss.Ny);
  W = zeros(1,ss.Nz);

  for idx=(1:length(shearZhat)-1)
    ss.shearZ_chat(idx) = shearZhat(idx) + abs(shearZhat(idx+1) - shearZhat(idx))/2;
    W(idx) = abs(shearZhat(idx) - shearZhat(idx+1));
  end
  for idx=(1:length(shearYhat)-1)
    ss.shearY_chat(idx) = shearYhat(idx) + abs(shearYhat(idx+1) - shearYhat(idx))/2;
    L(idx) = abs(shearYhat(idx) - shearYhat(idx+1));
  end

  L(end) = L(1);
  W(end) = abs(shearZhat(end-1) - shearZhat(end));

  % grid
  % non-hat terms are NxN
  shearZhat(end)=[]; shearYhat(end)=[]; shearXhat(end)=[];
  [shearY shearZ] = ndgrid(shearYhat, shearZhat);

  [shearY_c shearZ_c] = ndgrid(ss.shearY_chat, ss.shearZ_chat);

  shearX = zeros(size(shearY));
  shearX_c=zeros(size(shearY_c));

  % ensure necessary fields make it into `ss`
  ss.faultZ = faultZ;
  ss.shearY_c = shearY_c;
  ss.shearZ_c = shearZ_c;
  ss.shearY_chat = ss.shearY_chat;
  ss.shearZ_chat = ss.shearZ_chat;
  ss.shearZhat = shearZhat;

  disp('mesh done.')
  disp('beginning kernels...')

  % KERNELS
  % ---       General KVF PARAMETERS      ---
  c.command = 'compress';
  c.lambdaZ = ss.lambdaZ;
  c.tol = 1e-8;
  c.G = 30e3;
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % ---       s12 kernel stress on fault, by fault      ---
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/BP1v_ff-s12';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.Y = [faultX; faultY; faultZ'];
  c.X = [faultX_c; faultY_c; faultZ_c];
  c.L = Lf';
  c.W = Wf';

  kvf('Write', c.kvf, c, 4);
  disp('run these in a shell:')
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ff12 = c.write_hmat_filename;

  % ---       shear1212 kernel for shear-shear interaction      ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1212';
  c.write_hd_filename = './tmp/BP1v_ss-shear1212-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.Ny = ss.Ny;
  c.Nz = ss.Nz;

  c.transition = ss.transition;
  c.Y = [shearXhat'; ss.shearY_chat'; shearZhat'];
  c.X = [shearX_c(:)'; shearY_c(:)'; shearZ_c(:)'];

  c.Z = [shearX(:)'; shearY(:)'; shearZ(:)']; % for sizing purposes only
  % for some kernels, I found it easier to pass in the *_hat arrays rather than
  % the longer flattened arrays such as shearX_c(:)


  c.L = L;
  c.W = W;

  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;

  % ---       shear 1213 kernel for shear-shear interaction     ---
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1213';
  c.write_hd_filename = './tmp/BP1v_ss-shear1213-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear kernels
  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % ---       shear 1312 kernel for shear-shear interaction       ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1312';
  c.write_hd_filename = '/tmp/BP1v_ss-shear1312-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % ---       shear 1313 kernel for shear-shear interaction       ---
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1313';
  c.write_hd_filename = './tmp/BP1v_ss-shear1313-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1313 = c.write_hmat_filename;

  % ---       shear 1212 kernels for fault-shear interaction      ---
  c.X = [faultX_c; faultY_c; faultZ_c];
  c.Y = [shearXhat'; ss.shearY_chat'; shearZhat'];
  c.Z = [shearX(:)'; shearY(:)'; shearZ(:)']; % hm sizing purposes only

  c.L = L;
  c.W = W;

  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/BP1v_fs-shear1212';
  c.write_hd_filename = './tmp/BP1v_fs-shear1212-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1212 = c.write_hmat_filename;

  % ---       shear1312 kernels for fault-shear interaction     ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/BP1v_fs-shear1312';
  c.write_hd_filename = './tmp/BP1v_fs-shear1312-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1312 = c.write_hmat_filename;

  % ---         s12 kernel for shear-fault interaction      ---
  c.X = [shearX_c(:)'; shearY_c(:)'; shearZ_c(:)'];
  c.Y = [faultX; faultY; faultZ'];
  c.L = Lf';
  c.W = Wf';

  c.Bfro = 1e-8;

  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/BP1v_sf-s12';
  c.write_hd_filename = './tmp/BP1v_sf-s12-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf12 = c.write_hmat_filename;

  % ---       s13 kernel for shear-fault interaction      ---
  c.greens_fn = 'okadaS13';
  c.write_hmat_filename = './tmp/BP1v_sf-s13';
  c.write_hd_filename = './tmp/BP1v_sf-s13-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf13 = c.write_hmat_filename;

  disp('kernels done.')

  r.ss = ss;

end

% run simulation
function out = run(b)
  addpaths();

  ss = b.ss;

  % these are QDBIM style defns
  % shear wave speed (m/s)
  Vs = 3464;
  % density (kg/m^3)
  rho = 2670;
  % shear modulus (MPa)
  %G = rho*Vs^2/1e6;
  G = 30e3; ss.G = G;

  %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                                                      %
  %         F R I C T I O N   P A R A M E T E R S        %
  %                                                      %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

  % - QDBIM style

  % reference friction coefficient
  ss.fo=0.6*ones(size(ss.faultZ));

  % Dieterich-Ruina R+S frictional parameters (velocity-weakening friction)
  ss.a=1e-2+Ramp((ss.faultZ-15e3)/3e3)*(0.025-0.01);
  ss.b=0.015*ones(size(ss.faultZ));

  % effective normal stress (MPa)
  ss.sigma=50.0*ones(size(ss.faultZ));

  % characteristic weakening distance (m)
  ss.Drs=8e-3*ones(size(ss.faultZ));

  % plate rate (m/s)
  ss.Vpl=1e-9*ones(size(ss.faultZ));
  ss.Vpl_scalar = 1e-9;

  % reference slip rate (m/s)
  ss.Vo=1e-6*ones(size(ss.faultZ));

  % Radiation damping coefficient
  ss.eta = G./(2*Vs);

  % Estimates of some key parameters
  VWp = find(ss.a < ss.b); % VW region
  % Critical nucleation size ( h* = pi/2 G b D_rs / (b-a)^2 / sigma )
  hstar=min(pi/2*G*ss.Drs(VWp).*ss.b(VWp)./(ss.b(VWp)-ss.a(VWp)).^2./ss.sigma(VWp));

  % Quasi-static cohesive zone ( coh0 = 9/32 G D_rs /(b*sigma) )
  % Note that for this QD simulation the cohesive zone will not change,
  % which would not be the case for a fully dynamic simulation
  coh = min(9/32*pi*G*ss.Drs(VWp)./ss.b(VWp)./ss.sigma(VWp));

  % Estimate of recurrence time ( T ~ 5(b-a)*sigma / G * R/Vpl )
  Ti = 5*mean((ss.b(VWp)-ss.a(VWp)).*ss.sigma(VWp)).*0.5.*(ss.faultZ(VWp(end))- ...
    ss.faultZ(VWp(1)))./(G*mean(ss.Vpl(VWp)));

  % Print information about discretization
  fprintf('\nGrid size = %.2f (m)\n', ss.dz);
  fprintf('VW zone = %.2f (km)\n', (ss.faultZ(VWp(end))-ss.faultZ(VWp(1)))/1e3);
  fprintf('Critical nucleation size = %.2f (m)\n',hstar);
  fprintf('QS Cohesive zone = %.2f (m)\n',coh);
  fprintf('Est. Recurrence time = %.2f (yr)\n\n', Ti/3.15e7);

  %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                                                      %
  %                   R H E O L O G Y                    %
  %                                                      %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
  % Values taken for wet olivine - Hirth, G. and D. Kohlstedt (2003)

  % Driving strain rate (1/s)

  % solve for using the epsilon dote
  % equations 6 and 7 in the overleaf
  % add factor of plate rate out front
  Dv = ss.probW-ss.transition;
  Df = ss.lambdaZ;
  w = 10;

  % Power-Law Exponent
  ss.n = 3.0*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);
  n_scalar = 3.0;

  ss.e12p_plate = zeros(ss.Nz*ss.Ny, 1);
  ss.e13p_plate = zeros(ss.Nz*ss.Ny, 1);
  for i=1:1:ss.Ny
    for j = 1:1:ss.Nz
      x2 = ss.shearY_chat(i);
      x3 = ss.shearZ_chat(j);
      x2p = x2/(Dv);
      x3p = (x3-Df)/(Dv);

      summ12 = 0;
      summ13 = 0;
      m = 1;

      % also consider term as percentage of the total
      sterm12 = e12Terms(x2p, x3p, m, n_scalar, w);
      sterm13 = e13Terms(x2p, x3p, m, n_scalar, w);
      while sterm12 >= summ12*0.005 || sterm13 >= summ13*0.005;
        summ12 = summ12 + sterm12;
        summ13 = summ13 + sterm13;
        m = m + 1;
        sterm12 = e12Terms(x2p, x3p, m, n_scalar, w);
        sterm13 = e13Terms(x2p, x3p, m, n_scalar, w);
      end
      ss.e12p_plate((j-1)*ss.Ny+i) = ss.Vpl_scalar * ( 1/(2*w) + (1/w) * summ12);
      ss.e13p_plate((j-1)*ss.Ny+i) = ss.Vpl_scalar * ( (-1/(w*(n_scalar^0.5))) * summ13 );

    end
  end

  % flatten
  ss.e12p_plate = ss.e12p_plate';
  ss.e13p_plate = ss.e13p_plate';
  %ss.e12p_plate = ss.e12p_plate(:);
  %ss.e13p_plate = ss.e13p_plate(:);

  %ss.e12p_plate = 1e-14*ones(length(ss.shearY_chat)*length(ss.shearZ_chat),1);
  %ss.e13p_plate =      zeros(length(ss.shearY_chat)*length(ss.shearZ_chat),1);

  % Strengh profile

  ss.A = 10e-19;

  s120 = (ss.e12p_plate./ss.A).^(1./ss.n);
  s130 = (ss.e13p_plate./ss.A).^(1./ss.n);
  e120 = zeros(size(s120));
  e130 = zeros(size(s120));

  %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                                                       %
  %         N U M E R I C A L   S O L U T I O N           %
  %                                                       %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

  % state parameters
  ss.dgfF=4;
  ss.dgfS=4;
  %% Initialize State Vector
  Y0=zeros(ss.M*ss.dgfF+length(ss.shearY_chat)*length(ss.shearZ_chat)*ss.dgfS,1);

  % Fault patches
  % state vector is (slip; tau; log(theta Vo / D_rs); log(V / Vo) )
  Y0(1:ss.dgfF:ss.M*ss.dgfF) = zeros(ss.M,1);
  Y0(2:ss.dgfF:ss.M*ss.dgfF) = max(ss.a).*ss.sigma.*asinh(ss.Vpl./ss.Vo/2.* ...
    exp((ss.fo+ss.b.*log(ss.Vo./ss.Vpl))./max(ss.a))) + ss.eta.*ss.Vpl;
  Y0(3:ss.dgfF:ss.M*ss.dgfF) = ss.a./ss.b.*log(2*ss.Vo./ss.Vpl.* ...
    sinh((Y0(2:ss.dgfF:ss.M*ss.dgfF)-ss.eta.*ss.Vpl)./ss.a./ss.sigma))-ss.fo./ss.b;
  Y0(4:ss.dgfF:ss.M*ss.dgfF) = log(ss.Vpl./ss.Vo);

  % Shear zones
  Y0(ss.M*ss.dgfF+1:ss.dgfS:end)=s120;
  Y0(ss.M*ss.dgfF+2:ss.dgfS:end)=s130;
  Y0(ss.M*ss.dgfF+3:ss.dgfS:end)=e120;
  Y0(ss.M*ss.dgfF+4:ss.dgfS:end)=e130;

  % load kernels using hmmvp
  hm.ff12 = hmmvp('init', b.ff12, 4);
  hm.ss1212 = hmmvp('init', b.ss1212, 32);
  hm.ss1312 = hmmvp('init', b.ss1312, 32);
  hm.ss1213 = hmmvp('init', b.ss1213, 32);
  hm.ss1313 = hmmvp('init', b.ss1313, 32);
  hm.fs1212 = hmmvp('init', b.fs1212, 4);
  hm.fs1312 = hmmvp('init', b.fs1312, 4);
  hm.sf12 = hmmvp('init', b.sf12, 4);
  hm.sf13 = hmmvp('init', b.sf13, 4);

  % initialize the function handle with
  % set constitutive parameters
  yp=@(t,y) odeBP1v_hm(t,y,ss, hm);
  disp('begin solving...')
  tic
  % Solve the system
  options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6);
  [t,Y]=ode45_2(yp,[0 1*3.15e7],Y0,options);
  disp('done solving.')
  toc
  %%

  Y = Y'; % Y has slip and shear stress

  % Compute the instantaneous derivative
  Yp=zeros(length(t)-1,size(Y,2));
  for k=1:length(t)-1
      Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
  end

  %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                                                      %
  %                    F I G U R E S                     %
  %                                                      %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

  % save entire Y vector
  mat2np(Y, 'pickles/BP1vHM2_Y.pkl', 'float64');

  % Strain rate at center
  Ep=sqrt(Yp(:,ss.M*ss.dgfF+floor(length(ss.shearY_chat)/2)*ss.dgfS+3:ss.dgfS* ...
            length(ss.shearY_chat):end)'.^2 +...
          Yp(:,ss.M*ss.dgfF+floor(length(ss.shearY_chat)/2)*ss.dgfS+4:ss.dgfS* ...
            length(ss.shearY_chat):end)'.^2);
  mat2np(Ep, 'pickles/BP1vHM2_strainCenter.pkl', 'float64');


  % strain rate over whole ductile area
  Epall = sqrt( Yp(:,ss.M*ss.dgfF+3:ss.dgfS:end)'.^2 +...
               Yp(:,ss.M*ss.dgfF+4:ss.dgfS:end)'.^2);
  mat2np(Epall, 'pickles/BP1vHM2_strainAll.pkl', 'float64');

  % Velocity
  y.V = Yp(:,1:ss.dgfF:ss.M*ss.dgfF); % Slip rate (m/s)
  y.V = y.V';
  y.tau = Y(:,2:ss.dgfF:ss.M*ss.dgfF);            % Shear stress (MPa)
  mat2np(y.V, 'pickles/BP1vHM2_V.pkl', 'float64');

  % fault slip figure
  clf;
  imagesc(log10(y.V)); colorbar;
  title('Slip Rate')
  xlabel('time steps')
  ylabel('fault mesh block')
  saveas(gcf, 'figures/BP1vHM2_slip.png')

  % strain rate at center of ductile region
  clf;
  imagesc(log10(Ep)); colorbar;
  title('Strain Rate of Center of Ductile Region')
  xlabel('Time Steps')
  ylabel('Block')
  saveas(gcf, 'figures/BP1vHM2_strainCenter.png')

  % ---         Movies        ---
  Smovie=false;
  if Smovie
    disp('begin shear movie')
    clf;
    fig = figure;
    fname = 'figures/BP1vHM2_strain.gif';
    for idx = 1:size(Epall, 2)
      oneE = Epall(:,idx);
      oneEsq = reshape(oneE, [ss.Ny, ss.Nz]);
      imagesc(oneEsq'); colormap(hot); colorbar; caxis([])
      title(idx)
      drawnow
      frame = getframe(fig);
      im{idx} = frame2im(frame);

      [A,map] = rgb2ind(im{idx},256);
      if idx==1
        imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',0.1);
      else
        imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0.1);
      end

    end
  disp('shear movie done')
  end

  Fmovie=false;
  if Fmovie
    % velocity movie
    disp('begin fault movie')
    clf;
    fig = figure;
    fname='figures/BP1vHM2_slip.gif';
    for idx = 1:size(y.V,1)
      oneV = y.V(idx,:)';
      imagesc(oneV); colorbar; caxis([1e-15 5e-14]);
      title(idx)
      drawnow
      frame = getframe(fig);
      im{idx} = frame2im(frame);

      [A, map] = rgb2ind(im{idx}, 256);
      if idx==1
        imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',0.1);
      else
        imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0.1);
      end
    end
  disp('fault movie done')
  end

  out.hm = hm;

end

% XXX --------  Helper Functions ------------ XXX

% add paths to hmmvp
function addpaths()
  addpath('include/hmmvp/matlab')
  addpath('ODEsolving')
  addpath('include/')
end

% boxcar function
function bc = BC(x)
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
  bc = boxc(x);
end

% ramp function
function rmp = Ramp(x)
  rmp = x.*BC(x-1/2)+HS(x-1);
end

function h = HS(x)
  h = 0+x>=0;
end

% terms in the summation for e12
function y = e12Terms(x2p, x3p, m, n, w)
  t1 = cosh( (m*pi*(1-x3p))/(w*(n^0.5)) );
  t2 = cos( (m*pi*x2p)/(w*(n^0.5)) );
  t3 = cosh( (m*pi)/(w*(n^0.5)) );

  y = (t1/t3) * t2;
end

% terms in the summation for e13
function y = e13Terms(x2p, x3p, m, n, w)
  t1 = sinh( (m*pi*(1-x3p))/(w*(n^0.5)) );
  t2 = sin( (m*pi*x2p)/(w*(n^0.5)) );
  t3 = cosh( (m*pi)/(w*(n^0.5)) );

  y = (t1/t3) * t2;
end
