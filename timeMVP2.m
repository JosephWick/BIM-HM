% timeMVP.m

function varargout = tMVP (varargin)

    [varargout{1:nargout}] = feval(varargin{:});

end

% handles KVFs for all values of N, tolerances 1e-6 and 1e-8
function r = build()
  addpaths();

  Ns = [10, 20, 30, 40, 50];
  % total elems in mesh: [125, 1000, 15625, 125000, 1000000]
  probDim = 200e3; %200km
  transition = 40e3;

  disp('Run the following commands in a shell to build all HMs')

  i = 1;
  N = Ns(i);
  r.N = N;

  eps = 1e-12;
  nc = (-N/2:N/2);

  Xhat = linspace(0, probDim, N);
  Yhat = tan(nc*pi/(2.5*max(nc)))*32e3;
  Zhat = transition+tan((0:N)'*pi/(2.2*(N+eps)))*transition;

  L1 = abs(Xhat(1)-Xhat(2))*ones(1,N); % xhat direction is uniform meshing
  L2 = zeros(1,N);
  L3 = zeros(1,N);

  Xchat = Xhat+(L1(1)/2);
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

  length(Xhat);
  length(Yhat);
  length(Zhat);

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
  c.Y = [Xc; Yc; Zc];

  [L3f, L2f, L1f] = ndgrid(L3, L2, L1);
  L1f = L1f(:)';
  L2f = L2f(:)';
  L3f = L3f(:)';
  c.L = [L1f; L2f; L3f];

  %disp(size(c.X))
  %disp(size(c.Y))
  %disp(size(c.L))

  %  create kernel based on mesh
  G = 30e3;
  nu = 0.25; % this is fine
  % theta is zero
  % eps12 can be 1; others zero

  % general things
  c.command = 'compress';
  c.greens_fn = 'time';
  c.err_method = 'mrem-fro';
  c.allow_overwrite = 1;

  %c.Bfro = 1e-8;

  % tolerance 1e-6;
  c.tol = 1e-6;
  c.write_hmat_filename = './tmp/t3d_e6';
  c.write_hd_filename = './tmp/t3d_e6-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.kvfe6 = c.write_hmat_filename;

  % tolerance 1e-8;
  c.tol = 1e-8;
  c.write_hmat_filename = './tmp/t3d_e8';
  c.write_hd_filename = './tmp/t3d_e8-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.kvfe8 = c.write_hmat_filename;

end

function measure(b)
  N = b.N;
  disp(N)

  t_3d_6 = -1;
  t_3d_8 = -1;
  t_dense = -1;

  rng('default');
  v = zeros(N,N,N);
  ia = floor(N/4);
  ib = floor(3*N/4);
  v(ia:ib,ia:ib,ia:ib) = 1.0;
  v = v(:);

  numIter = 500;


  disp('timing e6...')
  t6 = [];
  for j=1:5
    hme6 = hmmvp('init', b.kvfe6);
    now = tic();
    x = hmmvp('mvp', hme6, v);
    t = toc(now);
    t6(j) = t;
  end
  t_3d_6 = mean(t6);

  disp('timing e8...')
  t8 = [];
  for j=1:5
    hme8 = hmmvp('init', b.kvfe8);
    now = tic();
    x = hmmvp('mvp', hme8, v);
    t = toc(now);
    t8(j) = t;
  end
  t_3d_8 = mean(t8);

  disp('timing dense...')
  td = [];
  for j = 1:5
    dense = hmmvp('extract', hme8, 1:1:N*N*N, 1:1:N*N*N);
    now = tic();
    x = dense*v;
    t = toc(now);
    td(j) = t;
  end
  t_dense = mean(td);

disp('RESULTS:')
disp('dense: ')
t_dense
disp('e6')
t_3d_6
disp('e8')
t_3d_8

end

% XXX --- Private --- XXX
function addpaths()
  addpath('include/hmmvp/matlab')
  addpath('ODEsolving')
  addpath('include/')
  addpath('./tmp/')
end
