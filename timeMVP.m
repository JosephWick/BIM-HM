% timeMVP.m

function varargout = time (varargin)

    [vararout{1:nargout}] = feval(varargin{:});

end

% handles KVFs for all values of N, tolerances 1e-6 and 1e-8
function r = build()
  addpaths();

  Ns = [10, 50, 100, 500, 1000];
  probDim = 200e3; %200km
  transition = 40e3;

  r.kvfs_2d_6 = [];
  r.kvfs_2d_8 = [];
  r.kvfs_3d_6 = [];
  r.kvfs_3d_8 = [];
  r.N = Ns;

  disp('Run the following commands in a shell to build all HMs')

  for i=1:length(Ns)
    N = Ns(i);
    nstring = string(N);

    eps = 1e-12;
    nc = (-N/2:N/2);

    Xhat = linspace(0, probDim, N);
    Yhat = tan(nc*pi/(2.5*max(nc)))*32e3;
    Zhat = transition+tan((0:N)'*pi/(2.2*(N+eps)))*transition;

    L1 = abs(Xhat(1)-Xhat(2))*ones(1,N); % xhat direction is uniform meshing
    L2 = zeros(1,N);
    L3 = zeros(1,N);

    Xchat = Xhat+L1(1);
    Ychat = zeros(1,N);
    Zchat = zeros(1,N);
    for idx=(1:length(Zhat)-1)
      L2(idx) = abs(Yhat(idx) - Yhat(idx+1));
      L3(idx) = abs(Zhat(idx) - Zhat(idx+1));

      Ychat(idx) = Yhat(idx) - abs(Yhat(idx) - Yhat(idx+1))/2;
      Zchat(idx) = Zhat(idx) - abs(Zhat(idx) - Zhat(idx+1))/2;
    end
    L2(end) = L2(1);
    L3(end) = abs(Zhat(end-1) - Zhat(end));
    Yhat(end)=[]; Zhat(end)=[];

    % --- 3D Mesh ---
    [Z, Y, X] = ndgrid(Zhat, Yhat, Xhat);
    X = X(:)';
    Y = Y(:)';
    Z = Z(:)';
    c.X = [X; Y; Z];

    [Zc,Yc,Xc] = ndgrid(Zchat, Ychat, Xchat);
    Xc = Xc(:)';
    Yc = Yc(:)';
    Zc = Zc(:)';

    [L3f, L2f, L1f] = ndgrid(L3, L2, L1);
    L1f = L1f(:)';
    L2f = L2f(:)';
    L3f = L3f(:)';

    %  create kernel based on mesh
    G = 30e3;
    nu = 0.25; % this is fine
    % theta is zero
    % eps12 can be 1; others zero

    disp('mesh done, making kernel')

    % src is upper left, rec is center
    kernel = zeros((N+1)^3, (N+1)^3);
    for i = 1:N^3
      for j = 1:N^3
        kernel(i,j) = computeStressVerticalShearZone_s12(...
        Xc(i), Yc(i), Zc(i), ...
        X(j), Y(j), Z(j), L1f(j), L2f(j), L3f(j), 0, ...
        0, 1, 0, 0, 0, 0, ...
        G, nu);
      end
    end

    disp('kernel done')

    % general things
    c.command = 'compress';
    c.greens_fn = 'time';
    c.err_method = 'mrem-fro';
    c.allow_overwrite = 1;

    csvwrite('K.csv', K);
    K = real(kernel);
    K(isnan(K)) = 0;
    c.K = K;

    c.Bfro = 1e-8;

    % tolerance 1e-6;
    c.tol = 1e-6;
    c.write_hmat_filename = './tmp/timingHM';
    c.write_hd_filename = './tmp/timingHM-hd';
    c.kvf = [c.write_hmat_filename '.kvf'];
    c.kvf
    kvf('Write', c.kvf, c, 4);
    cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
    disp(cmd)

    aaaaa

    % tolerance 1e-8;
    c.tol = 1e-8;
    c.write_hmat_filename = ['./tmp/t3d_n' nstring '_e8'];
    c.write_hd_filename = [c.write_hmat_filename '-hd'];
    c.kvf = [c.write_hmat_filename '.kvf'];
    kvf('Write', c.kvf, c, 32);
    cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
    disp(cmd)
    %r.kvfs_2d_8(length(r.kvfs_3d_8)+1) = c.write_hmat_filename;

  end

  b.Ns = Ns;

end

function measure(b)

  Ns = b.Ns;
  t_3d_6 = [];
  t_3d_8 = [];
  t_dense = [];

  Ns = [10, 25, 50, 100, 500, 1000];
  for i=1:length(Ns)
    N = Ns(i);
    v = ones(N*N*N)

    tees = [];
    for j=1:3
      hme6 = hmmvp('load', './tmp/timing_3d_n' + nstring + '_e6');
      now = tic();
      x = hmmvp('mvp', hme6, v);
      t = toc(now)
      tees(j) = t;
    end
    t_3d_6(i) = mean(tees);

    tees = [];
    for j=1:3
      hme8 = hmmvp('load', './tmp/timing_3d_n' + nstring + '_e8');
      now = tic();
      x = hmmvp('mvp', hme8, v);
      t = toc(now)
      tees(i) = t;
    end
    t_3d_8(i) = mean(tees);

    tees = [];
    for j = 1:3
      dense = hmmvp('extract', hme8, 1:1:N*N*N, 1:1:N*N*N);
      now = tic();
      x = dense*v;
      t = toc(now)
      tees(i) = t;
    end
    t_dense(i) = mean(tees);

  end

  T = [t_dense; t_3d_6; t_3d_8]

end

% XXX --- Private --- XXX
function addpaths()
  addpath('include/hmmvp/matlab')
  addpath('ODEsolving')
  addpath('include/')
  addpath('./tmp/')
end
