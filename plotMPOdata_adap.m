function plotMPOdata_adap(source,tmin,tlim,coloresid,cases, input, eigS )

switch coloresid
    case 1
        colores = ['k']
    case 2
        colores = ['k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r']; % 2
    case 3        
        colores = ['k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b'];  % 3
    case 4
        colores = [ [0,0,0],'r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m'];  % 4
       %colores = ['k','r','k','r','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m'];  % 4
    case 5
        colores = ['k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g']; % 5
    case 6
       colores = ['k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c']; % 5        
    case 8 
       colores = ['k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c']; % 5        
   case 9
       colores = ['k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c']; % 5        
        
end
                
     fig1=figure(1);clf(1);figure(1)
     fig2=figure(2);clf(2);figure(2)
     fig3=figure(3);clf(3);figure(3)
     fig4=figure(4);clf(4);figure(4)
    %fig5=figure(5);clf(5);figure(5)
     fig6=figure(6);clf(6);figure(6)
     fig7=figure(7);clf(7);figure(7)     
    %fig8=figure(8);clf(8);figure(8)     
     fig9=figure(9);clf(9);figure(9)
    %fig10=figure(10);clf(10);figure(10)
                        
    for qq=1:numel(source)
                
           %Aid = source{qq}.Aid;
           N =  source{qq}.N;
           dtfs = source{qq}.dtfs;            
           sysID = source{qq}.sysID;  
           vibrations_info = source{qq}.vibrations_info;  
           getreducedosc = source{qq}.getreducedosc;
           rho_out = source{qq}.rho_out;
           rhoscs = source{qq}.rhoscs;
           AccTruncerr = source{qq}.AccTruncerr;
           Truncerr = source{qq}.Truncerr;
           oide = source{qq}.oide;        
           Nb = source{qq}.Nb;
           Q = numel(Nb); 
           Nt     = source{qq}.Nt;
           chivals= source{qq}.chivals;            
           CSsites = diag(sysID);              
           %S = source{qq}.S;                                      
           %secstep = source{qq}.secstep;     
           select_oscs = source{qq}.getreducedosc; % assume every vib state has been calculated
           %eigS = source{qq}.eigS;  
           %repeats = source{qq}.repeats;                           
           %losses = source{qq}.losses;           
           %select_oscs = [1 2]
           %select_oscs =  [1:(Q-1)/repeats, Q];                                 
           %select_oscs = [1 5 10 15 20 ]
           
            tmaxfs = size(rho_out,1) * dtfs;
            tmax   = tmaxfs * 0.188/1000; % 0.188 = 1 ps       
            dt     = 0.188/1000 * dtfs;  
            nmax   = tmax/dt;         
            nmax   = floor(nmax);
            time = (0:nmax-1)*dt/0.188*1000;
            pop_idx=diag(sysID);                                  

            % Trace                            
                tr=zeros(nmax,1);                  
                for tt=1:nmax
                   for ii=1:N          
                      tr(tt)=tr(tt)+(rho_out(tt,pop_idx(ii)));                                
                   end
                end
                
            % Normalize
                rho_out=normalize_rho(rho_out);
                                            
            % Distance (after normalization)
                distance=zeros(nmax,1);                
                for tt=1:nmax
                   for ii=1:N          
                      distance(tt)=distance(tt)+(ii-1)* rho_out(tt,pop_idx(ii));                                
                   end
                end
                
%             % Integrated bond dimension
%             totchi = cell(numel(Aid),1);
%             maxchi = 1;
%             minchi = 1;
%             for i = 1:numel(Aid)                 
%                   totchi{i} = zeros(Nt,1);
%                   for tt=1:Nt
%                      %totchi{i}(tt) = sum( chivals{i}(tt,:),2 ); 
%                      totchi{i}(tt) = mean( chivals{i}(tt,:),2 ); 
%                   end
%                   maxchi = max( maxchi, max( totchi{i} ));
%                   minchi = min( minchi, min( totchi{i} ));
%             end
    
            if 1  
            % Plot Populations and Coherences
            figure(1)
                          
                % Coherences
                subplot(2,3,1)    
                plot(time,real(rho_out(:,sysID(~tril(ones(N))))),colores(qq), time, imag(rho_out(:,sysID(~tril(ones(N))))), colores(qq)); hold on    % plot only upper part of        xlabel('t','interpreter','latex')
                ylabel('$\rho_{ij}$','interpreter','latex')
                title('Coherences','interpreter','latex')
                xlim([0,tlim])
                ylim([-0.3,0.3])

                % Populations
                subplot(2,3,4)    
                plot(time,rho_out(:,diag(sysID)),colores(qq));hold on               
               %plot(time,rho_out(:,pop_idx([1 2 3 8 9])),colores(qq));hold on               
               %plot(time,rho_out(:,pop_idx([1 2 3])));hold on               
               %plot(time,rho_out(:,pop_idx([8 9 10])));hold on               
                
                title('Populations','interpreter','latex')
                xlabel('Time(fs)','interpreter','latex')
                ylabel('$\rho_{ii}$','interpreter','latex','Clipping','off')
                xlim([0,tlim])
                ylim([0,1])

                % distance, <x>
                subplot(2,3,2)            
                plot(time,distance',colores(qq));hold on                    
                title('$\chi$','interpreter','latex')
                xlabel('Time(fs)','interpreter','latex')
                title('$\langle \hat x \rangle$(t)','interpreter','latex')
                xlim([0,tlim])
                %ylim([0 5])
                %ylim([0 2])
                %ylim([0 12])
                grid on

                % trace
                subplot(2,3,5)            
                plot(time,tr,colores(qq)) ;hold on       
                xlabel('Time(fs)','interpreter','latex')
                title('${\rm{Tr}}{\hat{\rho}_S}$','interpreter','latex')
                xlim([0,tlim])
                %ylim([0,1.2])

                % accumulated error
                %subplot(2,3,3)                   
                %plot(time,AccTruncerr,colores(qq));hold on   
                %xlabel('$t[\rm{fs}]$','interpreter','latex') 
                %title('Accumulated Error','interpreter','latex')
                %xlim([0,tmaxfs])
                
                    % Populations at selected sites
                    subplot(2,3,3)                       
                    %plot(time, sum( rho_out(:,pop_idx([1 2])),2 ) ,'b','LineStyle','--','LineWidth',2);hold on
                    %plot(time, sum( rho_out(:,pop_idx(3)),2 ) ,'b','LineStyle','-','LineWidth',2);hold on
                    
                    %plot(time, sum( rho_out(:,CSsites(3:end)),2 ) ,colores(qq));hold on        
                    %plot(time, sum( rho_out(:,CSsites(4:end)),2 ) ,colores(qq),'LineStyle',':','LineWidth',2);hold on
                    
                     plot(time, sum( rho_out(:,CSsites(3:end)),2 ) ,colores(qq),'LineStyle','-','LineWidth',2);hold on
                     plot(time, sum( rho_out(:,pop_idx([1])),2),colores(qq),'LineStyle','--','LineWidth',1);
                     plot(time, sum( rho_out(:,pop_idx([2])),2),colores(qq),'LineStyle','--','LineWidth',1);
                     
                    
                     grid on
                    title('$\rho_{nn}$','interpreter','latex')
                    xlabel('Time(fs)','interpreter','latex')
                    xlim([tmin,tlim])

                % local compression error
                subplot(2,3,6)           
                plot(time,Truncerr,colores(qq));hold on
                xlabel('Time(fs)','interpreter','latex')
                title('Local Error','interpreter','latex')   
                xlim([0,tlim])  
            end

            
            
%             % Population profile at selected times
%             figure(10);
%             subplot(5,1,1)
%             stem( linspace(1,10,10), rho_out(50, pop_idx) ), hold on, title('Time=50 fs','interpreter','latex')
%             ylim([0 1])
%             grid on
%             subplot(5,1,2)
%             stem( linspace(1,10,10), rho_out(100, pop_idx) ), hold on, title('Time=100 fs','interpreter','latex')
%             ylim([0 1])
%             grid on
%             subplot(5,1,3)
%             stem( linspace(1,10,10), rho_out(200, pop_idx) ), hold on, title('Time=200 fs','interpreter','latex')
%             ylim([0 1])
%             grid on
%             subplot(5,1,4)
%             stem( linspace(1,10,10), rho_out(400, pop_idx) ), hold on, title('Time=400 fs','interpreter','latex')
%             ylim([0 1])
%             grid on
%             subplot(5,1,5)
%             stem( linspace(1,10,10), rho_out(1000, pop_idx) ), hold on, title('Time=50 fs','interpreter','latex')
%             ylim([0 1])
%             grid on
            
                
%             % Dynamics Contour         
%                 rho_exciton = zeros(size(rho_out,1) , size(rho_out,2) );
%                 for tt=1:nmax % First        
%                     % change basis
%                    rho_exciton(tt,:) = reshape( S'* reshape(rho_out(tt,:),N,N) *S , 1,[]);
%                 end
%                 
%                
%              % Population of delocalised excitons
%              figure(2)   
% 
%                  subplot(2,1,1)                 
%                 %plot(time,  sum( rho_exciton( :, pop_idx([2 3 4 5 6 8]) ), 2),colores(qq),'LineStyle','-');hold on               
%                  plot(time,  sum( rho_exciton( :, pop_idx([2 3 4 5 6 7 9 10]) ), 2),colores(qq),'LineStyle','-');hold on               
%                  title('Eigenstate populations','interpreter','latex')
%                  xlabel('Time(fs)','interpreter','latex')
%                  ylabel('$\rho_{ii}$','interpreter','latex','Clipping','off')
%                  xlim([0,tlim])
%                  %ylim([0,1])
%                
%                  subplot(2,1,2)
%                  plot(time, sum( abs(rho_out(:,sysID(~tril(ones(N))))), 2)  ,colores(qq)); hold on    
%                  ylabel('$\rho_{ij}$','interpreter','latex')
%                  title('$\sum |\rho_{nm}|, m < n$','interpreter','latex')
%                  xlim([0,tlim])
%                  %ylim([0,0.])
%              
             % Contour plot of population dynamics
             figure(3)
                subplot(121)
                levels =  [0:0.025:0.4];
                ro=abs(rho_out(:,diag(sysID)));
                [C,h]=contourf(linspace(1,N,N),time,ro, levels,'ShowText', 'on' );                
                set(h,'LineColor','none')
                
                    %imagesc( linspace(1,N,N),time,ro, [0 0.4] )
                    %pcolor( linspace(1,N,N),time,ro )
                    %shading( 'flat' )
                
                xlabel('Site Index (n)','interpreter','latex')
                ylabel('Time (fs)','interpreter','latex')
                title('Site basis','interpreter','latex')
                caxis([0 0.4])
                colorbar()
                ylim([0 tlim])
                xlim([0.5 N+0.5])                                                

                subplot(122)
                levels =  [0:0.025:0.4];
                %ro=abs( rho_exciton(:,diag(sysID)) );                
                %[C,h]=contourf(eigS,time,ro, levels, 'ShowText', 'on' );                
                %set(h,'LineColor','none')
                
                    %imagesc( linspace(1,N,N),time,ro, [0 0.4] )                
                    %uimagesc( diag(eigS),time,ro, [0 0.4] )                
                    %uimagesc( eigS,time,ro, [0 0.4] )                
                    %pcolor( eigS,time,ro )
                    %shading( 'flat' )
                
                xlabel('Energy (cm$^{-1})$','interpreter','latex')
                ylabel('Time (fs)','interpreter','latex') 
                title('Exciton basis','interpreter','latex')
                colorbar()
                %caxis([0 0.4])
                ylim([0 tlim])
                %xlim([1 N])
    
           % losses
%            figure(10)
%                                 
%                 plotAid= N+1; % column corresponding to pop of site 2
% 
%                 subplot(211)
%                     plot(time, losses(:,1), colores(qq) ), hold on
%                     plot(time, losses(:,1).*rho_out(:,1), colores(qq) ,'LineStyle','--'), hold on
%                     grid on
%                     title('$O_{11}$','interpreter','latex')   
%                     xlim([0,tlim])  
% 
%                 subplot(212)
%                     plot(time, losses(:, plotAid ), colores(qq) ), hold on
%                     plot(time, losses(:, plotAid ).*rho_out(:, source{qq}.sysID(2,2) ), colores(qq) ,'LineStyle','--'), hold on
%                     grid on
%                     title('$O_{22}$','interpreter','latex')   
%                     xlim([0,tlim])  
%                    % ylim([-0.1,1])
         

                
%            % Plot coherence between first and 8,9th site    
%            figure(4)                      
%                 subplot(211)
%                 plot(time,imag(rho_out(:,8)),colores(qq) ), hold on 
%                 xlim([0,tlim])       
%                 ylim([-0.1 0.1])                         
%                 xlabel('t','interpreter','latex')
%                 ylabel('$\rho_{ij}$','interpreter','latex')
%                 title('$\rho_{18}$','interpreter','latex')
%                 
%                 subplot(212)
%                 plot(time,imag(rho_out(:,9)),colores(qq), 'LineStyle','-'), hold on             
%                 %plot(time,imag(rho_out(:,10*2 + 2)),colores(qq), 'LineStyle',':'), hold on             
%                 %plot(time,abs(rho_out(:,10*2 + 1)),colores(qq), 'LineStyle','-.'), hold on                             
%                 %plot(time,real(rho_out(:,2)),colores(qq) ,time, imag(rho_out(:,2)),colores(qq)), hold on 
%                 %plot(time,real(rho_out(:,3)),colores(qq) ,time, imag(rho_out(:,3)),colores(qq),'LineStyle','--'), hold on                             
%                 %ylim([-0.2,0.3])    
%                 %ylim([0 0.3])       
%                 xlim([0,tlim])       
%                 %ylim([-0.1 0.1])                         
%                 xlabel('t','interpreter','latex')
%                 ylabel('$\rho_{ij}$','interpreter','latex')
%                 title('$\rho_{19}$','interpreter','latex')

           
                
                %             figure(2)            
                % 
                %                 % Dipole autocrrelation function
                %                 subplot(1,3,1)
                %                 %corr = sum(rho_out(:,2:N),2);
                %                 corr = rho_out(:,4) + rho_out(:,7);                
                % 
                %                  plot(time,real(corr),colores(qq) ,time, -imag(corr),colores(qq)), hold on 
                %                 %plot(input(:,1)* 0.05308837 * 1000, input(:,3),'k--', input(:,1)* 0.05308837 * 1000, input(:,4),'k--');      
                %                 xlim([0,tlim])
                %                 ylim([-1,2])    
                %                 xlabel('t','interpreter','latex')
                %                 ylabel('$\rho_{ij}$','interpreter','latex')
                %                 title('Correlation Function','interpreter','latex')
                % 
                %                 % Coherences
                %                 subplot(1,3,2)
                %                 plot(time,real(rho_out(:,8))/max(real(rho_out(:,8))),colores(qq) ,time, imag(rho_out(:,8))/max(real(rho_out(:,8))),colores(qq)), hold on 
                %                 %plot(input(:,1)* 0.05308837 * 1000, input(:,9)/max(input(:,9)),'k--', input(:,1)* 0.05308837 * 1000, input(:,10)/max(input(:,9)),'k--');    
                %                 xlim([0,tlim])       
                %                 ylim([-1,1])    
                %                 xlabel('t','interpreter','latex')
                %                 ylabel('$\rho_{ij}$','interpreter','latex')
                %                 title('$\rho_{12}$','interpreter','latex')
                % 
                % 
                %                 % Populations
                %                 subplot(1,3,3)
                %                 plot(time,real(rho_out(:,5)) * 2/3,colores(qq) ,time, real(rho_out(:,9))* 2/3,colores(qq)), hold on 
                %                 %plot(input(:,1)* 0.05308837 * 1000, input(:,6),'k--', input(:,1)* 0.05308837 * 1000, input(:,7),'k--'); 
                %                 xlim([0,tlim])
                %                 ylim([-1,1])    
                %                 xlabel('t','interpreter','latex')
                %                 ylabel('$\rho_{ij}$','interpreter','latex')
                %                 title('$\rho_{12}$','interpreter','latex')

             
            % Vibrational observables (all states)    
            if 0    
            if vibrations_info    
                % Plot dynamics of <n|rho|n>
               figure(5);clf(5);figure(5)
            else
               clf(5);figure(5)
            end                                         
                %clf(3);figure(3)
                orderplot=zeros(Q,N);
                cc = 1;
                for ss=1:Q
                for nn=1:N                            
                    orderplot(ss,nn) = cc;
                    cc = cc+1;
                end
                end
                for oo=1:numel(select_oscs)

                    whichosc=oide(getreducedosc(oo));
                    n_idx = diag( reshape(find(ones(Nb(whichosc),Nb(whichosc))),[Nb(whichosc),Nb(whichosc)]) );        
                    rhovtrace = sum( rhoscs{select_oscs(oo)}(:,n_idx) ,2); 
                    for it=1:Nt
                        rhoscs{select_oscs(oo)}(it,:) = rhoscs{select_oscs(oo)}(it,:) / sum(rhoscs{select_oscs(oo)}(it,n_idx),2);       
                        %rhoscs{select_oscs(oo)}(it,:) = rhoscs{select_oscs(oo)}(it,:) %/ sum(rhoscs{select_oscs(oo)}(it,n_idx),2);       
                    end                  
                    %orderplot = reshape(reshape(find(ones(Q,N)),Q,N)',[],1);
                    subplot(Q,N,orderplot(oo))

                    plot(time,rhoscs{select_oscs(oo)}(:,n_idx), colores(qq), 'LineWidth',1.5),hold on
                    plot(time, rhovtrace ,'k--', 'LineWidth',2);
                    grid on
                    xlabel('Time(fs)','interpreter','latex')
                    ylabel('$<n|\hat{\rho}|n>$','interpreter','latex')
                    ylim([-0.01 0.1])
                    xlim([0 tlim])
                    title(sprintf('osc. %i', oo),'interpreter','latex')                    
                end
            end
            
            % Plot osc_i-osc_j oscillators of a given site
            if 0       
            if vibrations_info    
                                               
                % Pseudomodes (frequency first)
                whichsite  = 2;               
                listofoscs = [1:N:N*Q] + (whichsite-1); % oscillators belonging to site whichsite
                
                % Surr. (site first)
                %whichsite   = 1;               
                %listofoscs = [ (whichsite-1)*Q ; % oscillators belonging to site whichsite
                                
                omin=1;
                omax=10; % max Q
                rangeosc = listofoscs(omin:omax);                
                
                % Plot dynamics of <n|rho|n>
                figure(6);                                                                            
                for oo=1:numel( rangeosc )
                    
                    whichosc=oide(rangeosc(oo));                     
                    n_idx= diag( reshape(find(ones( Nb(whichosc),Nb(whichosc)) ),[ Nb(whichosc), Nb(whichosc) ]) );        
                    rhovtrace = sum( rhoscs{rangeosc(oo)}(:,n_idx) ,2); 
                    for it=1:Nt % normalize
                        rhoscs{rangeosc(oo)}(it,:) = rhoscs{rangeosc(oo)}(it,:) / sum(rhoscs{rangeosc(oo)}(it,n_idx),2);       
                    end
                
                    subplot( numel(rangeosc) ,1, oo )

                    plot(time,rhoscs{ rangeosc(oo) }(:,n_idx), colores(qq), 'LineWidth',1.5),hold on
                    plot(time, rhovtrace ,'k--', 'LineWidth',2);
                    grid on
                    xlabel('Time(fs)','interpreter','latex')
                    ylabel('$<n|\hat{\rho}|n>$','interpreter','latex')
                    %ylim([0 0.4])
                    xlim([0 tlim])
                    title(sprintf('osc. %i', rangeosc(oo)),'interpreter','latex')                        
                end
                
            end    
            end
            
            % Vibrational observables (only selected ones), Assuming select_osc contains all oscs and repeats
            if 0
            vibengy = zeros(Nt, N);    
            if vibrations_info    
                
                % Plot dynamics of <n|rho|n>
                figure(5);
                
                jij = 1;
                 
                hfidx = [Q:Q:N*Q];
                %hfidx(1) = 1;
                
                for oo=1:numel(hfidx) 
                                        
                        whichosc = oide( hfidx(oo) );

                        n_idx= diag( reshape(find(ones(Nb(whichosc),Nb(whichosc))),[Nb(whichosc),Nb(whichosc)]) );        
                        % normalize
                        rhovtrace = sum( rhoscs{ hfidx(oo) }(:,n_idx) ,2); 
                        for it=1:Nt
                            rhoscs{ hfidx(oo) }(it,:) = rhoscs{hfidx(oo)}(it,:) / sum(rhoscs{hfidx(oo)}(it,n_idx),2);       
                        end
                        % Vibrational energy                                                                                               
                        vibengy(:,jij) =  abs( rhoscs{ hfidx(oo) }(:,n_idx)) * [0:1210:(Nb(whichosc)-1)*1210]';                                                  
                                                                                  
                        subplot(N,1,jij)

                           %plot(time,rhoscs{  hfidx(oo) }(:,n_idx), colores(qq),'LineWidth',1),hold on
                           %plot(time,rhoscs{  hfidx(oo) }(:,n_idx),'LineWidth',1),hold on                            
                            plot(time,rhoscs{  hfidx(oo) }(:,n_idx([1 2 3 4 5])),'LineWidth',2),hold on
                            
                            plot(time, rhovtrace ,colores(qq), 'LineWidth',2,'LineStyle','--');
                            grid on
                            xlabel('Time(fs)','interpreter','latex')
                            ylabel('$<n|\hat{\rho}|n>$','interpreter','latex')
                            ylim([0 1])
                            xlim([0 tlim])
                            title(sprintf('osc. %i', oo),'interpreter','latex')
                            ylim([0 0.4])   
                        
                            jij = jij + 1;                                                                                        
                end
                   
                % Vibrational energy at each site
                figure(6)                                                                                                              
                        subplot(1,2,1)
                       %levels = 0 : 5 : max(vibengy(1:tlim,:));
                       %[C,h]=contourf( linspace(1,N,N), time, vibengy, levels, 'ShowText', 'on');                                       
                       %set(h,'LineColor','none')
                        imagesc(  linspace(1,N,N), time, vibengy ) 
                        colormap( hot )
                        ylabel('Time (fs)','interpreter','latex')                             
                        colorbar()
                        ylim([ 0 tlim])
                         
                        subplot(1,2,2)
                        %levels = 0 : 5 : 350;
                        %[C,h]=contourf( linspace(1,N,N), time, vibengy, levels, 'ShowText', 'on');                
                        %[C,h]=contourf( linspace(1,N,N), time, vibengy, 'ShowText', 'on');                                        
                        %set(h,'LineColor','none')
                        %caxis([0 500])
                        imagesc(  linspace(1,N,N), time, vibengy ) 
                        ylabel('Time (fs)','interpreter','latex')                             
                        colorbar()                        
                        ylim([ 0 tlim])    
                        caxis([0 150])
                        
                figure(7)
                    subplot(211)
                        %plot(time, sum(vibengy,2) / max( sum(vibengy,2)), colores(qq) ), hold on
                        %plot(time, distance / max(distance), colores(qq),'LineStyle',':' ), hold on
                        %plot(time, -abs( rho_out(:,2) / max(rho_out(:,2)) ), colores(qq) ), hold on
                        plot(time, vibengy(:,[2 7 8 9]), colores(qq) ), hold on
                        
                        
                        ylabel('$\langle \hat H_v \rangle$','interpreter','latex')
                        xlabel('Time (fs)','interpreter','latex')                                                     
                    subplot(212)
                        plot(time, sum(vibengy,2), colores(qq) ), hold on
                        xlabel('Time (fs)','interpreter','latex')          
                        
            end
            end
                                                                   
            % Bond dimensionS, O_11 and O_NN
            if 0
                figure(8)
                ax(1)=subplot(211)                
                    [C,h]=contourf(linspace(1,N*Q-1,N*Q-1),time, chivals{1},  'ShowText', 'on' );                
                    %h.LevelList=linspace(1,20,20);   
                    ylim([0 tlim])                                
                    ylabel('Time (fs)','interpreter','latex')                                
                    xticks( linspace(1,N*Q-1,N*Q-1) )
                    xticklabels( mod( linspace(1,N*Q-1,N*Q-1), Q))                
                    caxis([1 35])
                    colorbar()


                    ax(2)=subplot(212)
                    [C,h]=contourf(linspace(1,N*Q-1,N*Q-1),time, chivals{end},  'ShowText', 'on' );                                                                                    
                    %h.LevelList=linspace(1,20,20);  
                    ylim([0 tlim])
                    ylabel('Time (fs)','interpreter','latex')
                    xticks( linspace(1,N*Q-1,N*Q-1) )
                    xticklabels( mod( linspace(1,N*Q-1,N*Q-1), Q))                          
                    caxis([1 35])
                    colorbar()

                    bonds=[ mod( linspace(1,N*Q-1,N*Q-1),4); linspace(1,N*Q-1,N*Q-1)]
            end
            
            % Bond dimensions (integrated)
            if 0        
            figure(8);            
            for i = 1:numel(Aid)
                [m,n]=find( sysID==Aid(i) );                
                subplot(N,N, (m-1)*N+n )
                plot( time, totchi{i} , colores(qq), 'LineWidth', 1), hold on
                grid on
                xlim([0 tlim])
                ylim([minchi maxchi])                
            end 
            end
            
            if 0
            % Bond dimensions (upper diagonal)
            figure(8);
            disp('preparing bond dimension dynamics plots')
            for i = 1:numel(Aid)

                [m,n]=find( sysID==Aid(i) );                
                chivals{i}(chivals{i}==0)=1;
                chivals{i}(isnan(chivals{i}));

                subplot(N,N, (m-1)*N+n )
                
                    bonds = repmat( linspace(1,Q/repeats,Q/repeats), 1, repeats * N); bonds(end)=[]; 

                    title('$\chi(t)$','interpreter','latex')       
                    [C,h]=contourf( linspace(1,N*Q-1,N*Q-1), time, chivals{i} / chivals{i}(1,1) ); hold on                        
                    %surfc(linspace(1,N*Q-1,N*Q-1),time,chivals{i})
                    %xticklabels(bonds)
                    %surfc(linspace(1,N*Q-1,N*Q-1),time,chivals{i})
                    %h.LevelList=linspace(0,0.6,20);                       
                    ylabel('Time (fs)','interpreter','latex')
                    caxis([1 20])
                    colorbar()
                    ylim([0 tlim])                                        

            end             
            [ mod( linspace(1,N*Q-1,N*Q-1),4); linspace(1,N*Q-1,N*Q-1)]
             end
            
            if 0
            % Bond dimensions (single column)
            figure(8);
            for i = 1:numel(Aid)

                [m,n]=find( sysID==Aid(i) );   

                chivals{i}(chivals{i}==0)=1;

                subplot(numel(Aid),1,i)

                title('$\chi(t)$','interpreter','latex')       
                [C,h]=contourf(linspace(1,N*Q-1,N*Q-1),time, chivals{i} / chivals{i}(1,1) );                        
                %surfc(linspace(1,N*Q-1,N*Q-1),time,chivals{i})
                %h.LevelList=linspace(0,0.6,20);                       
                ylabel('Time (fs)','interpreter','latex')
                colorbar()
                ylim([0 tlim])
            end 
            end
            
           
            %            % Fourier Transform (in cm-1) (MAKE SURE rho_out HAS BEEN NORMALISED!)
            %                 %input = rho_out(:,2:N); % |g><1| + |g><2| + ... 
            %                 ftinput = rho_out(:,4)+rho_out(:,7);
            %                 ftinput(isnan(ftinput))=0;
            %         
            %                 resolution=2000;
            %                 omega=linspace(-700,700,resolution);
            %                 timecm=(0:nmax-1)*dt;
            %                 Aw = zeros(1,resolution);
            %                 for w=1:resolution
            %                     %for it=1:size(ftinput(~isnan(ftinput)),1)
            %                      for it=1:nmax
            %                         %Aw(w) = Aw(w) + omega(w)* exp(+1j*omega(w)*timecm(it)) * (ftinput(it));
            %                         Aw(w) = Aw(w) + exp(-1j*omega(w)*timecm(it)) * ftinput(it);           
            %                    end
            %                 end
            %                 disp('Fourier fisnihed')
            %         
            %             % Absorption plot (in eV)
            %             figure(4);
            %         
            %                 plot(omega,real(Aw)/max(real(Aw)),colores(qq));hold on
            %                 %plot((input(:,2)*100 -12405.05 )*2*pi,input(:,5)/max(input(:,5)), 'k--');hold on
            %                 xlim([-700,700])
            %                 grid on      
            %                         
            %                 end
 
%             % Simulation time
%             figure(9);%              
%                 plot(time, secstep, colores(qq)), hold on
%                 xlim([0,tlim])
%                 title('Simulation Time','interpreter','latex')
%                 xlabel('Time (fs)','interpreter','latex')
%                 grid on
                 
            if qq ~= numel(source)                                  
                
                    disp('Click Me Baby One more Time')
                    pause         
            else
                            
     
                
    end                        
        disp('Finished plotting data')
end