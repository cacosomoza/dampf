
% % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                   %
% DAMPF, Alejandro Somoza and Oliver Marty, (2019)  %  
% Instute of Theoretical Physics, Ulm Universität   %
%                                                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % %

clear, disp('clearing, start')

% Some physical constansts
unit_cm_fs = 0.00018836515661233488;  % fs to cm-1
eVtocm=8065.6;                        % 1 eV = 8065.6 cm-1
kB = 8.6173303e-5 * eVtocm;           % cm.-1 / K

% Solver Parameters
tmaxfs = 200;                         % in fs
dtfs = 1.0;                           % timestep (in fs)

tmax = tmaxfs*unit_cm_fs;             % fs/1000 = ps * 0.188 cm-1 / 1ps
dt   = unit_cm_fs * dtfs;             % numerical timestep (cm-1)Â´
Nt = round(tmax/dt);                  % number of total steps

% Global parameters
parallel=0;

savedata=0;                           % save output in a .mat file
epsilon=1e-16;                        % Elements of the superoperator smaller than epsilon will be discarded
normalization= 0;                     % norm 1, 0
vibrations_info = 1;                  % Calculate vibrational observables
displayinforate = 1;                  % display information every 
holdplots   = 1;                      % superpose plots            

% Bond dimension envelope   
comp_tol = 1 - 1e-8;
chi0 = 1;
chimax = 6;

% Model parameters
N  = 3;        % number of sites
J  = 500;      % Exciton Coupling
for i=1:N                     
Omega(i)=0;    % site energies            
end
T  = 0;        % Temperature, cm-1 (it will be set to 0 if surrogate oscillators are used)

% Dephasing
deph   = 0;     % dephasing ON/OFF
gdphfs = 60;    % dephasing rate (fs)

 
% Initial state System
rhoS = zeros(N,N);
rhoS(1,1) = 1; % local OK

%environment_type = 'pseudomodes';
 environment_type = 'surrogate';

 inputcase = '2_osc';

 switch inputcase
    
    case '1_osc'
        ws = [1500];                            % Osc. frequency
        Q = numel(ws);                          % number of oscillators per site
        Nb = [4];                               % truncation number
        HR = [0.20]*1.0;                        % Huang Rhys factor
        gdmpfs = [100];                         % osc. damping rates  (fs)
        damp = [1];                             % damping yes or no
        g = ones(1,Q-1) * 0.0;                  % phonon-phonon coupling

        T  = 200;                               % in wavenumbers
        T  = T .* ones(1,Q);           
        Q = numel(ws);                          % number of oscillators per site    
        getreducedosc = linspace(1,Q*N,Q*N);    % get reduced state of oscillators indicated in the list    
        
        for qq=1:Q
            jc(qq) = sqrt(HR(qq))*ws(qq);                     % Exciton-Phonon coupling strength 
            gdmp(qq) =  damp(qq) * 1/unit_cm_fs/gdmpfs(qq);   % osc. damping rates  (1/cm-1)
        end
         
    case '2_osc'
        ws = [1500,500];                        % Osc. frequency
        Q = numel(ws);                           % number of oscillators per site
        Nb = [4,4];                              % truncation number
        HR = [0.20 0.10]*1.0;                    % Huang Rhys factor
        gdmpfs = [100 50];                       % osc. damping rates  (fs)
        damp = [1 1];                            % damping yes or no
        g = ones(1,Q-1) * 100.0;                   % phonon-phonon coupling

        T  = 0;
        T  = T .* ones(1,Q);           
        Q = numel(ws);                           % number of oscillators per site    
        getreducedosc = linspace(1,Q*N,Q*N);     % get reduced state of oscillators indicated in the list    
        
        for qq=1:Q
            jc(qq) = sqrt(HR(qq))*ws(qq);                     % Exciton-Phonon coupling strength 
            gdmp(qq) =  damp(qq) * 1/unit_cm_fs/gdmpfs(qq);   % osc. damping rates  (1/cm-1)
        end
        
    case '3_osc'
        ws = [1500,1000,500];                    % Osc. frequency
        Q = numel(ws);                           % number of oscillators per site
        Nb = [3,4,5];                            % truncation number
        HR = [0.10 0.10 0.10]*1.0;               % Huang Rhys factor
        gdmpfs = [100 100 20];                    % osc. damping rates  (fs)
        damp = [1 1 1];                          % damping yes or no
        g = ones(1,Q-1) * 250.0;                  % phonon-phonon coupling

        T  = 0;
        T  = T .* ones(1,Q);    
        Q = numel(ws);                           % number of oscillators per site    
        getreducedosc = linspace(1,Q*N,Q*N);     % get reduced state of oscillators indicated in the list    
        
        for qq=1:Q
            jc(qq) = sqrt(HR(qq))*ws(qq);                     % Exciton-Phonon coupling strength 
            gdmp(qq) =  damp(qq) * 1/unit_cm_fs/gdmpfs(qq);   % osc. damping rates  (1/cm-1)
        end    
end

disp( inputcase );
             
% Some auxiliary variables labeling properties of ALL oscillators (N*Q)
oide=[]; % gives the 1:Q index from all the oscillators 1:N*Q
we=[];
Nbe=[];
gdmpe=[];

% If w_{site,q} The ordering of the oscillators seems to perform (error-wise) better with an ordering where equal q's lie together.
% However, when performing calculations with surrogate osc., it is necessary to use a local ordering with oscs. from the same site together.

switch environment_type 
    case 'pseudomodes'
        together = 0; 
    case 'surrogate'
        together = 1;
end

if together
    % w1 w2 w1 w2    
    oide = repmat( linspace(1,Q,Q), 1,N );
    we   = repmat( ws, 1, N);
    Nbe  = repmat( Nb, 1, N);    
    gdmpe= repmat(gdmp,1, N);    
else
    % % w1 w1 w2 w2
    for qq=1:Q    
        oide  = [oide qq*ones(1,N)];        % for every osc. index labels which site belongs to
        we    = [we  ws(qq)*ones(1,N)];
        Nbe   = [Nbe Nb(qq)*ones(1,N)];        
        gdmpe = [gdmpe gdmp(qq)*ones(1,N)];
    end
end

oide
we
Nbe
gdmpe

% System Hamiltonian: Tight binding (single-excitation manifold) 
gdph  = deph *1/unit_cm_fs/gdphfs *ones(1,N);   % spin dephasing rate (1/cm-1) 
HS   = zeros(N, N);

jarr = J*ones(N-1, 1);
HS   = diag(Omega) + (diag(jarr, 1) + diag(jarr, -1));
HS2  = HS*HS;
[S,eigS]=eig(HS);

% System propagator with dephasing
Id = eye(N);
H = kron( eye(N), HS ) - kron( HS.', eye(N) ); 
D = zeros(N^2, N^2);
for j = 1 : N
     pe= Id(:,j)*Id(j,:);    
     D =  D + gdph(j) * ( kron( pe,pe ) - 0.5*kron( Id,pe ) - 0.5*kron(pe,Id) );
end
Us = reshape( expm( dt * (-1i*H + D) ), N, N, N, N );

% System indeces
sysID=reshape(find(ones(N,N)),[N,N]);                        % 1d idx of rho_S in matrix form
transmap=[reshape(sysID,[N^2 1]) reshape(sysID.',[N^2 1])];  % Side by side 1d idx of sysID and its transpose
loweridx=find(tril(ones(N,N)-eye(N,N)));                     % indexes of lower-diagonal elements
pop_idx=diag(sysID);   
i=0;

for m=1:N
for n=m:N
    i=i+1;
    Aid(i) = sub2ind([N N],m,n); % one dimensional index of |m><n|
end
end
TERMS = numel(Aid);

% % Reservoir occupation numberf
meanresT = @(x) 1/( exp(x) - 1 );
         
% Outputs
rho_out     = zeros(Nt,N^2); % reduced system state
Truncerr    = zeros(1,Nt);   % local compression error
AccTruncerr = zeros(1,Nt);   % acc. error
occup       = zeros(Nt,N*Q); % <N>   
occup2      = zeros(Nt,N*Q); % <N^2> 
Hb_aveg     = zeros(Nt,1);   % <Hv>

chivals = cell(numel(Aid),1);
for i=1:numel(Aid)
   chivals{i} = ones(Nt,N*Q-1); 
end

% Vibrational states
if vibrations_info
    rhoscs      = cell(numel(getreducedosc),1); 
    for oo=1:numel(getreducedosc)
        oscid = getreducedosc(oo);
        rhoscs{oo} = zeros(Nt,Nbe(getreducedosc(oo))^2); 
    end
end

%% Basic operators
disp('Basic Operators')

DD = cell(Q,1);
DDb= cell(Q,1);
a  = cell(Q,1);
ad = cell(Q,1);

maxdim=200;

a_safe  = sqrt(diag(1:maxdim-1, 1));
ad_safe = a_safe';
    
  for qq=1:Q    

    % Displacement Operator       
    DD{qq}  = expm( -1i * dt * (jc(qq) * a_safe + conj(jc(qq)) * a_safe')) ;
    % Free Bath Evolution
    DDb{qq} = expm( -1i * dt * ws(qq) * ad_safe * a_safe );    
    % Annhilation / Creation Operators
    a{qq}  =  a_safe;
    ad{qq} = ad_safe;
    
    % Truncate operators
    a{qq}   =   a{qq}( 1:Nb(qq),1:Nb(qq) );
    ad{qq}  =  ad{qq}( 1:Nb(qq),1:Nb(qq) );
    DD{qq}  =  DD{qq}( 1:Nb(qq),1:Nb(qq) );
    DDb{qq} = DDb{qq}( 1:Nb(qq),1:Nb(qq) );
    
  end
Id = eye(N);

%% Basis for Bath Operators

disp('Basis for Bath Operators')

% Fock Basis
deltaop = cell(Nb(qq)^2,1);
for i=1:Nb(qq)^2
    deltaop{i} = zeros(Nb(qq),Nb(qq));
    deltaop{i}(i) = 1;
end
    
% New Basis
basis=cell(Q,1);
sgntrnsp=cell(Q,1);
for qq=1:Q
           
    % identity (the rest are traceless)    
    basis{qq} = cell(Nb(qq)^2,1);
    basis{qq}{1} = speye(Nb(qq))/sqrt(Nb(qq));
    
    % diagonal
    for j = 1:Nb(qq)-1
        basis{qq}{j+1} = spalloc(Nb(qq),Nb(qq),Nb(qq));
        basis{qq}{j+1}(j,j) = 1/sqrt(2);
        basis{qq}{j+1}(j+1,j+1) = -1/sqrt(2);
    end
    
    % off-diagonal
    c = Nb(qq)+1;
    for j = 1:Nb(qq)
        for k = j+1:Nb(qq)
            basis{qq}{c} = spalloc(Nb(qq),Nb(qq),Nb(qq));
            basis{qq}{c}(j,k) = 1/sqrt(2);
            basis{qq}{c}(k,j) = 1/sqrt(2);
            c = c + 1;
            basis{qq}{c} = spalloc(Nb(qq),Nb(qq),Nb(qq));
            basis{qq}{c}(j,k) = 1/sqrt(2);
            basis{qq}{c}(k,j) = -1/sqrt(2);
            c  = c + 1;        
        end
    end
    
    % QR decomposition (Basis transformation) , this makes the operators orthogonal to each other
        % Vectorize {o} and store in BT
        BT = [];
        for j = 1:length(basis{qq})
            % reshape(input,[],n) automatically adjusts the number of rows
            BT = [BT,reshape(basis{qq}{j},[],1)]; 
        end
        % Apply QR
        [BT,~] = qr(BT);
        
        % When taking the complex adjoint, some op. spit out a -1, others a 1. (hermitian, antihermitian)
        sgntrnsp{qq} = zeros( Nb(qq)^2, 1);
        
        for j = 1:Nb(qq)^2
            % Reshape back into operators
            basis{qq}{j} = reshape(BT(:,j),Nb(qq),Nb(qq));                        
            sgntrnsp{qq}(j) = issymmetric(basis{qq}{j});            
        end
        
end

%% Basis for 2-body operators

basis2 = cell(Q-1,1);
for qq = 1:Q-1
    basis2{qq} = cell(Nb(qq)^2, Nb(qq+1)^2);
    for jj = 1:Nb(qq)^2
    for kk = 1:Nb(qq+1)^2
        basis2{qq}{jj,kk} = kron( basis{qq}{jj}, basis{qq+1}{kk} );
    end
    end
end
 
%% One-body operators

disp('One-body operators')

% Write operators (in matrix form) in the {o} basis
% trace (o_j^\dag A o_k) = j, k element of A with respect to the basis {o}   
Idosc      = cell(Q,1);
DDleft     = cell(Q,1);
DDrght     = cell(Q,1);
DDswch     = cell(Q,1);
BBath      = cell(Q,1);
Damp_down  = cell(Q,1);
Damp_up    = cell(Q,1);
NN         = cell(Q,1);
NN2        = cell(Q,1);

a_transf  = cell(Q,1);
ad_transf = cell(Q,1);

for qq = 1:Q
    
    Idosc{qq}      = speye(Nb(qq),Nb(qq)) ;
    DDleft{qq}     = sparse(Nb(qq),Nb(qq));    
    DDrght{qq}     = sparse(Nb(qq),Nb(qq));    
    DDswch{qq}     = sparse(Nb(qq),Nb(qq));    
    BBath{qq}      = sparse(Nb(qq),Nb(qq));    
    Damp_down{qq}  = sparse(Nb(qq),Nb(qq));    
    Damp_up{qq}    = sparse(Nb(qq),Nb(qq)); 
    NN{qq}         = sparse(Nb(qq),Nb(qq)); 
    NN2{qq}        = sparse(Nb(qq),Nb(qq));     
    a_transf{qq}   = sparse(Nb(qq),Nb(qq)); 
    ad_transf{qq}  = sparse(Nb(qq),Nb(qq)); 

    for j = 1:length(basis{qq}) % 1:Nb{qq}^2
    for k = 1:length(basis{qq}) % 1:Nb{qq}^2
                                                  
        %  Occupation number
            NN{qq}(j,k)  =  trace( basis{qq}{j}'* ad{qq}*a{qq} *basis{qq}{k} );   
            
        %  Occupation number^2
            NN2{qq}(j,k) =  trace( basis{qq}{j}'* ad{qq}*a{qq}*ad{qq}*a{qq} *basis{qq}{k} );   
            
        % Interaction term
            % Multiplying on the left 
            DDleft{qq}(j, k) = trace( basis{qq}{j}' *(DD{qq}*basis{qq}{k}) );            
            
            % Multiplying on the right
            DDrght{qq}(j, k) = trace( basis{qq}{j}' *(basis{qq}{k}*DD{qq}') );            
            
            % Sandwhich with j=k
            DDswch{qq}(j, k) = trace( basis{qq}{j}' *DD{qq}*basis{qq}{k}*DD{qq}' );            
            
        % Bath free evolution
            BBath{qq}(j, k)  = trace( basis{qq}{j}' *DDb{qq}*basis{qq}{k}*DDb{qq}' );            
            
        % damping    
            Damp_down{qq}(j,k) = trace( reshape( basis{qq}{j},[Nb(qq)^2,1] )' * expm( dt*gdmp(qq)* (meanresT(ws(qq)/T(qq))+1) *( kron(a{qq},a{qq})   - 0.5*kron(Idosc{qq},ad{qq}*a{qq}) - 0.5*kron(ad{qq}*a{qq},Idosc{qq}) ) ) * reshape(basis{qq}{k},[Nb(qq)^2,1]) );                              
            Damp_up{qq}(j,k)   = trace( reshape( basis{qq}{j},[Nb(qq)^2,1] )' * expm( dt*gdmp(qq)* (meanresT(ws(qq)/T(qq))  ) *( kron(ad{qq},ad{qq}) - 0.5*kron(Idosc{qq},a{qq}*ad{qq}) - 0.5*kron(a{qq}*ad{qq},Idosc{qq}) ) ) * reshape(basis{qq}{k},[Nb(qq)^2,1]) );                  
            
    end
    end
    
end

%% Two body operators

disp('Two-body operators')

U_osc_int = cell(Q-1,1);

for bb=1:Q-1 %labels bonds (both even and odd, every site has U and Ubar)
        
        ad_a{bb} = kron( ad{bb} ,  a{bb+1} );
        a_ad{bb} = kron( a{bb}  , ad{bb+1} );  
        
        Osc_int{bb} =  expm( -1i * dt * g(bb) * (ad_a{bb} + a_ad{bb}) );  % dimension Nb(bb) X Nb(bb+1)                                               
        Osc_int_adj{bb} = Osc_int{bb}';
        
        U_osc_int{bb} = zeros( Nb(bb)^2, Nb(bb+1)^2, Nb(bb)^2, Nb(bb+1)^2 );
        U_osc_int{bb} = zeros( Nb(bb)^2, Nb(bb+1)^2, Nb(bb)^2, Nb(bb+1)^2 );
        
        for j1 = 1:Nb(bb)^2
        for j2 = 1:Nb(bb+1)^2
            for i1 = 1:Nb(bb)^2
            for i2 = 1:Nb(bb+1)^2
                U_osc_int{bb}(j1,j2,i1,i2) = trace( basis2{bb}{j1,j2}' * Osc_int{bb} * basis2{bb}{i1,i2} * Osc_int_adj{bb} );
            end
            end
        end
        end                          
end

%% Writing initial state as MPO 

disp('Initial state')

% Oscillators
osc0 = cell(Q,1);
c0   = cell(Q+1,1);  
partitionfunction = zeros(Q,1);
         
for qq = 1:Q
           
    if T(qq) == 0
        % osc. gnd state
        osc0{qq} = zeros(Nb(qq),Nb(qq));        
        osc0{qq}(1,1) = 1;                
    else
        % Thermal state
        osc0{qq} = zeros(Nb(qq),Nb(qq));        
        partitionfunction(qq) = 0;        
        for nn=1:Nb(qq)
            partitionfunction(qq) = partitionfunction(qq) + exp( -(nn-1)*ws(qq)/T(qq) );            
        end
        for nn=1:Nb(qq)
            osc0{qq}(nn,nn) = exp( -(nn-1)*ws(qq)/T(qq) )/partitionfunction(qq);
        end    

    end               
    
    % coefficients of |g><g| or |beta><beta| in the o{} basis
    c0{qq}=zeros(1,Nb(qq)^2);    
    for j = 1:Nb(qq)^2 
        c0{qq}(j) = trace( osc0{qq}*basis{qq}{j} );        
    end

    disp('Initial State of Oscillators')    
    osc0{qq}

end

% MPOs state

  % A0:  N X N cells (MPOs). Each N cells (A-matrices). Each ( - X - X Nb^2) array
  % OBC: For s=1 (1 X chi X Nb^2), for s=N (chi X 1 X Nb^2) and for other s (chi X chi X Nb^2)
    A0 = cell(numel(Aid),1); 
    for j = 1:numel(Aid) % sys        
        for l = 1:N*Q % osc. id 
            for m = 1:Nbe(l)^2 % local osc. basis           

                % Open Boundary Conditions on the MPO
                % Returns 1 X chi for l=1; chi X 1 for l=N and chi X chi for other ls               
                A0{j}{l}(:,:,m) = zeros(chi0^(l~=1),chi0^(l~=N*Q));         

                % Initial condition: Only non-zero MPO is the j-k block at t0
                %A0{j, k}{l}(1,1,m) = (j==site0)*(k==site0)*c0(m);                               
                A0{j}{l}(1,1,m) = rhoS(Aid(j))^(1/(N*Q)) * c0{oide(l)}(m);                           
                
            end
        end        
    end  
    A = A0;

% Blank MPO
BlankMPO = cell(1,N);
for k = 1:N*Q 
    BlankMPO{k} = zeros(1,1,Nbe(k)^2);
end
 

tracefactor=1;
    for qq=1:Q
        tracefactor = tracefactor * trace(basis{qq}{1})^N;
    end
    
%% Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('MPO Dynamics')

tic

for it = 1:Nt
    
    % Output reduced density matrix       
    for m=1:numel(Aid)       
        % Tracing out over osc. dof only implies o_1 (the rest are traceless)
        % Multiply A-matrices for a given m,n MPO.                       
        exp_temp = tracefactor * A{m}{1}(:,:,1);
        for j = 2:N*Q 
            exp_temp = exp_temp  * A{m}{j}(:,:,1);
        end
        rho_out(it,Aid(m)) = exp_temp;               
    end
                
    % Normalization
    if normalization    
    normfactor=sum(rho_out(it,pop_idx));    
    for i=1:numel(Aid)            
        A{i}= MPO_scalar(1/normfactor, A{i});        
    end        
    end
    
    % trace
    trloc=0;
    for ii=1:N
        trloc = trloc+rho_out(it ,pop_idx(ii));
    end
       
    if vibrations_info
    rhoscs = trace_system( A, Aid, it, tracefactor, rhoscs,  N, Q, oide, Nb, Nbe, getreducedosc, pop_idx, basis, parallel );
    end % End of Vibrational state analysis
    
    
    % System-bath interaction:
    A = evolution_ev( A, together, Aid, N, Q, Nb, Nbe, DDswch, DDleft, DDrght, parallel );
               
    % Vibrational Evolution and Damping   
    A = freebath_damping( A ,together, Aid,  N, Q, damp, oide, BBath, Damp_down, Damp_up, parallel );
        
    % Vibrational interactions ( requires w1 w2 w1 w2 order ) (ONLY FOR TSO)
    [A, bonderr] = oscosc_coupling( A, N, Q, Aid, U_osc_int, Nb, Nbe, comp_tol, parallel ,g, chimax);
    
    % System                               
    [A, truncerr, chivals] = system_evolution( A, it, Aid, epsilon, N, Q, Us, oide, Nbe, sgntrnsp, comp_tol, chivals, parallel, chimax );
    Truncerr(it) = sqrt(2*truncerr) + sqrt(2*bonderr);
    AccTruncerr(it) = AccTruncerr((it-1)^(it > 1)) + Truncerr(it);        

                                            
    % fill lower diagonal
    rho_out(it,:) = reshape( reshape(rho_out(it,:),[N,N]) + conj( reshape(rho_out(it,sysID'),N,N) - diag(rho_out(it,pop_idx)) ) ,1,[] );
    
    vibtracecheck = 0;
    for oo=1:numel(getreducedosc)
        whichosc=oide(getreducedosc(oo));
        n_idx= diag( reshape(find(ones(Nb(whichosc),Nb(whichosc))),[Nb(whichosc),Nb(whichosc)]) );
        vibtracecheck = vibtracecheck +  sum(rhoscs{oo}(it,n_idx),2);
    end
    
     % Time Tracker
    if mod(round(it), displayinforate)==0
       fprintf(' %f sec/step, ',toc/it);fprintf( 'Simulation Time: %.1f fs, trace:%.4f, local error: %d, Acc. error: %d,vibracheck: %.3f \n',(it)*dtfs, trloc, Truncerr(it),AccTruncerr(it),vibtracecheck  )       
    end
    
    if mod(round(it), 100) == 0

        try 
            save('data_all','-v7.3')
        catch
            disp('something went wrong while saving data_all.mat')
        end

        try
	       save('data','rho_out' ,'Aid','tmaxfs','dtfs','sysID',...
                'AccTruncerr','Truncerr','oide','getreducedosc','Nb','rhoscs','occup2','occup','Nt','chivals','N','vibrations_info','repeats','secstep','HS','S')
        catch
            disp('something went wrong while saving data.mat');
        end
        
    end
    

end % end of dynamics

totaltime=toc

if savedata
save('data')
end

%% Plot Populations and Coherences
    
    tmax   = tmaxfs * 0.188/1000; % 0.188 = 1 ps
    dt     = 0.188/1000 * dtfs;  
    nmax   = tmax/dt;             
    nmax   = floor(nmax);
    time    = (0:nmax-1)*dt/0.188*1000;    
                    
                    
    tlim=tmaxfs;
    pop_idx=diag(sysID);                                          
    
    ruby=rho_out;
    
    % Distance
    distance=zeros(nmax,1);    
    tr=zeros(nmax,1);    
    for tt=1:nmax
       for ii=1:N          
          tr(tt)=tr(tt)+(ruby(tt,pop_idx(ii)));          
          distance(tt)=distance(tt)+(ii-1)*(ruby(tt,pop_idx(ii)));          
       end
    end
    
    rho_out=normalize_rho(rho_out);
    
% Compact plot 
figure(1)
if holdplots
       
else
   clf(1);figure(1)
end
    % Coherences
    subplot(2,3,1)    
     plot(time,real(rho_out(:,sysID(~tril(ones(N))))), time, imag(rho_out(:,sysID(~tril(ones(N)))))); hold on    % plot only upper part of rho
    xlabel('Time(fs)','interpreter','latex') 
    ylabel('$\rho_{ij}$','interpreter','latex')
    title('Coherences','interpreter','latex')
    xlim([0,tlim])
    ylim([-0.5,0.5])

    % Populations
    subplot(2,3,4)    
    plot(time,rho_out(:,diag(sysID)));hold on               
    title('Populations','interpreter','latex')
    xlabel('Time(fs)','interpreter','latex') 
    ylabel('$\rho_{ii}$','interpreter','latex','Clipping','off')
    xlim([0,tlim])
    ylim([0,1])

    subplot(2,3,2)              
    title('$\chi(t)$','interpreter','latex')       
    plot(time,distance); hold on                    
           
            
    subplot(2,3,5)            
    plot(time,tr,'k') ;hold on       
    xlabel('Time(fs)','interpreter','latex') 
    title('${\rm{Tr}}{\hat{\rho}_S}$','interpreter','latex')
    tix=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tix,'%.4f'))
    xlim([0,tlim])
    
    subplot(2,3,3)                   
    plot(time,AccTruncerr,'r');hold on   
    xlabel('Time(fs)','interpreter','latex') 
    title('Accumulated Error','interpreter','latex')
    xlim([0,tlim])
    
    subplot(2,3,6)           
    plot(time,Truncerr,'k');hold on
    xlabel('Time(fs)','interpreter','latex') 
    title('Local Error','interpreter','latex')   
    tix=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tix,'%.4f'))
    xlim([0,tlim])  
   
% Excitonic Basis

    figure(2);clf(2);figure(2)
       
    %HS = zeros(N, N); jarr = J*ones(N-1, 1); HS = eye(N)* 0 + (diag(jarr, 1) + diag(jarr, -1)); [S,eigS]=eig(HS);    
    
    rho_exciton = zeros(size(rho_out,1) , size(rho_out,2) );
    for tt=1:nmax % First        
        % change basis
       rho_exciton(tt,:) = reshape( S'*reshape(rho_out(tt,:),N,N)*S , 1,[]);
    end
      
     % Contour
     figure(2)
     
        subplot(121)
        ro=abs(rho_out(:,diag(sysID)));
        [C,h]=contourf(linspace(1,N,N),time,ro);                          
        ylabel('Time (fs)','interpreter','latex')
        colorbar()
        ylim([0 tlim])
        xlim([1 N])
        
        subplot(122)
        ro=abs(rho_exciton(:,diag(sysID)));
       [C,h]=contourf(linspace(1,N,N),time,ro);                      
        ylabel('Time (fs)','interpreter','latex')
        colorbar()
        ylim([0 tlim])
        xlim([1 N])
                
%      % Coherences
%     subplot(2,1,1)    
%     plot(time,real(rho_exciton(:,sysID(~tril(ones(N))))), time, imag(rho_out(:,sysID(~tril(ones(N)))))); hold on    % plot only upper part of rho
%     xlabel('t','interpreter','latex')
%     ylabel('$\rho_{ij}$','interpreter','latex')
%     title('Coherences','interpreter','latex')
%     xlim([0,tlim])
%     ylim([-0.5,0.5])
% 
%     % Populations
%     subplot(2,1,2)    
%     plot(time,rho_exciton(:,diag(sysID)));hold on               
%     title('Populations','interpreter','latex')
%     xlabel('t[fs]','interpreter','latex')
%     ylabel('$\rho_{ii}$','interpreter','latex','Clipping','off')
%     xlim([0,tlim])
%     ylim([0,1])
    
% Bond dimensions
figure(5);clf;figure(5)
          
    for i = 1:numel(Aid)

        [m,n]=find( sysID==Aid(i) );

        chivals{i}(chivals{i}==0)=1;

        subplot(N,N, (m-1)*N+n )
            levels =  [1:1:20];
            title('$\chi(t)$','interpreter','latex')       
            [C,h]=contourf(linspace(1,N*Q-1,N*Q-1),time, chivals{i} / chivals{i}(1,1), levels, 'ShowText', 'on'  ); hold on                        
            %surfc(linspace(1,N*Q-1,N*Q-1),time,chivals{i})
            %h.LevelList=linspace(0,20,30);          
            caxis([1 chimax])
            ylabel('Time (fs)','interpreter','latex')
            colorbar()
            ylim([0 tlim])                                        

    end 
    
    
if 0
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 5];
    print('-f5',strcat('feedback'),'-dpng','-r400');
    % print('-f1',strcat('trdistance'),'-dpng','-r400');  
end
      
% Vibrational observables
if vibrations_info
    
% Plot dynamics of <n|rho|n>

   figure(3);
   if holdplots
       
   else
       clf(3);figure(3)
   end

    for oo=1:numel(getreducedosc)

        whichosc=oide(getreducedosc(oo));
        n_idx= diag( reshape(find(ones(Nb(whichosc),Nb(whichosc))),[Nb(whichosc),Nb(whichosc)]) );        
        rhovtrace = sum( rhoscs{oo}(:,n_idx) ,2); 
        subplot(numel(getreducedosc),1,oo)        

            plot(time,rhoscs{oo}(:,n_idx),'LineWidth',2),hold on
            plot(time, rhovtrace ,'k--', 'LineWidth',2);
            grid on
            xlabel('Time(fs)','interpreter','latex') 
            ylabel('$<n|\hat{\rho}|n>$','interpreter','latex')
            title(sprintf('osc. %i', oo),'interpreter','latex')

    end
end

% save an image
if 0
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 4 3];
        %saveas(gcf,'mandel_parameter','epsc');
        print -depsc2 -tiff -r300 -painters mandel_parameter.eps
        print('-f6','mandel_parameter','-dpng','-r400');
end


