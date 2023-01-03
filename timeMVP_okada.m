% timeMVP_okada.m

function varargout = tMVP (varargin)

    [varargout{1:nargout}] = feval(varargin{:});

end

% build
function r = build()
  addpaths();

  % parameters we need for okada fortan code
  % alpha
  % x,y,z of observer
  % D source depth
  % dip angle in degrees (0)
  % POT1, POT2, POT3, POT4 potentials
  % UX, UY, UZ, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ

  probDim = 200e3; %200km
  transition = 40e3;

  disp('Run the following commands in a shell to build all HMs')

  i = 4;
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

  [L3f, L2f, L1f] = ndgrid(L3, L2, L1);
  L1f = L1f(:)';
  L2f = L2f(:)';
  L3f = L3f(:)';

  %  create kernel based on mesh
  G = 30e3;
  nu = 0.25; % this is fine
  % theta is zero
  % eps12 can be 1; others zero

  % general things
  c.command='compress';
  c.greens_fn = 'okada';

end

% XXX --- Private --- XXX
function addpaths()
  addpath('include/hmmvp/matlab')
  addpath('ODEsolving')
  addpath('include/')
  addpath('./tmp/')
end
