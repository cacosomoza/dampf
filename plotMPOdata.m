function plotMPOdata(source,tmin,tlim,coloresid,cases)

switch coloresid
    case 2
        colores = ['k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r','k','r']; % 2
    case 3        
        colores = ['k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b','k','r','b'];  % 3
    case 4
        colores = ['k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m','k','r','b','m'];  % 4
    case 5
        colores = ['k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g','k','r','b','m','g']; % 5
end
        
    %colores = ['r','b','g','b','r','g','b','r','g','b','r','g','b']
    %colores = ['k','r','r','r','b','b','b','g','g','g','m','m','m']  
    %colores = ['k','r','r','r','b','b','b','g','g','g','m','m','m']

    
     fig1=figure(1);clf(1);figure(1)
     fig2=figure(2);clf(2);figure(2)
     fig3=figure(3);clf(3);figure(3)
     %fig4=figure(4);clf(4);figure(4)
     fig5=figure(5);clf(5);figure(5)
     
     
     
    for qq=1:numel(source)
                
        fprintf('chi: %i',source{qq}.chivals(1))
        sprintf('Nb:')
        source{qq}.Nb'
        
        N =  source{qq}.N;
        dtfs = source{qq}.dtfs;
        tmaxfs = source{qq}.tmaxfs;
        sysID = source{qq}.sysID;
        CSsites = diag(sysID);
        %CSsites = CSsites(end-5:end);
        %CSsites = CSsites([3:10]);

        tmin=tmin;
        tlim=tlim;

        sprintf('case: %i',cases(qq))   
        
        rho_out = source{qq}.rho_out;
        AccTruncerr = source{qq}.AccTruncerr;
        Truncerr = source{qq}.Truncerr;
        %secstep = source{qq}.secstep;
        oide = source{qq}.oide;        
        Nb = source{qq}.Nb;
        Nt     = source{qq}.Nt;
        chivals= source{qq}.chivals;
                    
        tmax   = tmaxfs * 0.188/1000; % 0.188 = 1 ps
        dt     = 0.188/1000 * dtfs;  
        nmax   = tmax/dt;             % number of total steps
        nmax   = floor(nmax);
        time = (0:nmax-1)*dt/0.188*1000;

        pop_idx=diag(sysID);  % 1d idx of populations
        
        % trace (before normalization)
        tr = sum( rho_out(:,pop_idx), 2); % first calculate trace to show decay   
        
        % normalization 
        rho_out = normalize_rho(rho_out); % then normalize to calculate <X> and every other variable
        
        % Distance
        distance=zeros(nmax,1);                
        for tt=1:nmax
           for ii=1:N                    
              distance(tt)=distance(tt)+(ii-1)*(rho_out(tt,pop_idx(ii)));          
           end
        end
           
%          figure(1)
%         % Coherences
%         subplot(2,1,1)    
%          plot(time,real(rho_out(:,sysID(~tril(ones(N))))),colores(qq), time, imag(rho_out(:,sysID(~tril(ones(N))))),colores(qq)); hold on    % plot only upper part of rho
%         xlabel('t','interpreter','latex')
%         ylabel('$\rho_{ij}$','interpreter','latex')
%         title('Coherences','interpreter','latex')
%         xlim([0,tlim])
%         ylim([-0.5,0.5])
% 
%         % Populations
%         subplot(2,1,2)    
%         plot(time,rho_out(:,diag(sysID)),colores(qq));hold on                                       
%         title('Populations','interpreter','latex')
%         xlabel('t[fs]','interpreter','latex')
%         ylabel('$\rho_{ii}$','interpreter','latex','Clipping','off')
%         xlim([0,tlim])
%         ylim([0,1])
%         ylim([0,0.42])
% 
%         % local error
%         figure(2)
%         subplot(3,1,1)
%         plot(time,Truncerr,colores(qq));hold on
%         xlabel('$t[\rm{fs}]$','interpreter','latex') 
%         title('Local Error','interpreter','latex')   
%         xlim([0,tlim])
%                 
%         % Trace and 
%         subplot(3,1,3)
%         plot(time,tr,colores(qq)) ;hold on         
%         xlabel('$t[\rm{fs}]$','interpreter','latex') 
%         title('${\rm{Tr}}{\hat{\rho}_S}$','interpreter','latex')
%         xlim([0,tlim])
%         ylim([0 1])
%         
%         % Acccumulated error                
%         subplot(3,1,2)      
%         plot(time,AccTruncerr,colores(qq));hold on   
%         xlabel('$t[\rm{fs}]$','interpreter','latex') 
%         title('Accumulated error','interpreter','latex')
%         xlim([0,tlim])
%         %ylim([0 1])
        
        
    % Compact plot 
        figure(1)
        % Coherences
        subplot(2,3,1)    
         plot(time,real(rho_out(:,sysID(~tril(ones(N))))),colores(qq), time, imag(rho_out(:,sysID(~tril(ones(N))))),colores(qq)); hold on    % plot only upper part of rho
        xlabel('t','interpreter','latex')
        ylabel('$\rho_{ij}$','interpreter','latex')
        title('Coherences','interpreter','latex')
        xlim([tmin,tlim])
        ylim([-0.5,0.5])

        % Populations
        subplot(2,3,4)    
        plot(time,rho_out(:,diag(sysID)),colores(qq));hold on                                       
        title('Populations','interpreter','latex')
        xlabel('t[fs]','interpreter','latex')
        ylabel('$\rho_{ii}$','interpreter','latex','Clipping','off')
        xlim([tmin,tlim])
        ylim([0,1.0])       

        % Distance
        subplot(2,3,2)            
        plot(time,distance,colores(qq)); hold on           
        title('$<\hat x>$','interpreter','latex')
        xlabel('$t[\rm{fs}]$','interpreter','latex')     
        xlim([tmin,tlim])
        %ylim([0 1.2])
        grid on

        % Trace and Acccumulated error
        subplot(2,3,5)            
        plot(time,tr,colores(qq)) ;hold on 
        plot(time,AccTruncerr,colores(qq));hold on   
        xlabel('$t[\rm{fs}]$','interpreter','latex') 
        title('${\rm{Tr}}{\hat{\rho}_S}$ \& Accumulated error','interpreter','latex')
        xlim([tmin,tlim])
        ylim([0 1])

%         % Bond dimension
%         subplot(2,3,3)                       
%         plot(time,chivals,colores(qq));hold on       
%         title('$\chi(t)$','interpreter','latex')
%         xlabel('$t[\rm{fs}]$','interpreter','latex')     
%         xlim([0,tlim])
        
        
        % Populations at selected sites
        subplot(2,3,3)                       
        plot(time,sum(rho_out(:,CSsites),2),colores(qq));hold on        
        title('$\rho_{nn}$','interpreter','latex')
        xlabel('$t[\rm{fs}]$','interpreter','latex')     
        xlim([tmin,tlim])
        

        % local error
        subplot(2,3,6)
        plot(time,Truncerr,colores(qq));hold on
        xlabel('$t[\rm{fs}]$','interpreter','latex') 
        title('Local Error','interpreter','latex')   
        xlim([tmin,tlim]) 


        %figure(2)
        % 
        %    subplot(2,1,1)
        %    plot(time,secstep,colores(qq)); hold on    
        %    subplot(2,1,2)
        %    plot(time,ones(1,size(time,2)) * sum(secstep)/3600,colores(qq)); hold on


  
    % Error    
        figure(2)          
        % Trace and Acccumulated error
        subplot(2,2,1)            
        plot(time,tr,colores(qq)) ;hold on 
        plot(time,AccTruncerr,colores(qq));hold on   
        xlabel('$t[\rm{fs}]$','interpreter','latex') 
        title('${\rm{Tr}}{\hat{\rho}_S}$ \& Accumulated error','interpreter','latex')
        xlim([tmin,tlim])
        ylim([0 1])

        % local error
        subplot(2,2,2)        
        plot(time,Truncerr,colores(qq));hold on
        xlabel('$t[\rm{fs}]$','interpreter','latex') 
        title('Local Error','interpreter','latex')   
        xlim([tmin,tlim])  
        
        if qq > 1
            % local error
            subplot(2,2,3)            
            cla()
            plot(time, abs( real(source{qq}.rho_out(:,diag(sysID))) - real(source{qq-1}.rho_out(:,diag(sysID)))) );hold on
            xlabel('$t[\rm{fs}]$','interpreter','latex') 
            title('Population diff.','interpreter','latex')   
            xlim([tmin,tlim])
            ylim([0,0.05])
            
            % local error
            subplot(2,2,4)
            cla()
            plot(time, abs(real(source{qq}.rho_out(:,sysID(~tril(ones(N))))) - real(source{qq-1}.rho_out(:,sysID(~tril(ones(N))))) ));hold on            
            xlabel('$t[\rm{fs}]$','interpreter','latex') 
            title('Coherence diff.','interpreter','latex')   
            xlim([tmin,tlim]) 
            ylim([0,0.05])
        end

    % Population and coherence        
   figure(3)     
        %Coherences
        subplot(1,2,1)    
        plot(time,real(rho_out(:,sysID(~tril(ones(N))))),colores(qq), time, imag(rho_out(:,sysID(~tril(ones(N))))),colores(qq)); hold on    % plot only upper part of rho
        xlabel('Time(fs)','interpreter','latex')
        ylabel('$\rho_{ij}$','interpreter','latex')
        title('Coherences','interpreter','latex')
        xlim([tmin,tlim])
        ylim([-0.15,0.15])
        %ylim([-0.4,0.4])

        % Populations
        subplot(1,2,2)    
        plot(time,rho_out(:,diag(sysID)),colores(qq));hold on                                       
        title('Populations','interpreter','latex')
        xlabel('Time(fs)','interpreter','latex')
        ylabel('$\rho_{ii}$','interpreter','latex','Clipping','off')
        xlim([tmin,tlim])
        ylim([0.0,0.2])       

%         % Contour
%         figure(2)
%         subplot(1,numel(cases),qq)
% 
%             ro=abs(rho_out(:,diag(sysID)));
%             [C,h]=contourf(linspace(1,N,N),time,ro);
%             h.LevelList=linspace(0,0.6,20);                       
%             ylabel('fs')
%             colorbar()
%             ylim([0 tlim])
%             xlim([1 N])
            
    if source{qq}.vibrations_info
        
        if 1
            % Vibrational observables        
            getreducedosc = source{qq}.getreducedosc;
            rhoscs = source{qq}.rhoscs;
            occup2 = source{qq}.occup2;
            occup  = source{qq}.occup;

            avgN = zeros(Nt,numel(getreducedosc));
            avgN2 =zeros(Nt,numel(getreducedosc));

            figure(5);%clf(5);figure(5)
            
            %lastosc=numel(getreducedosc);
             lastosc=5;
            for oo=1:lastosc
                whichosc=oide(getreducedosc(oo));
                n_idx= diag( reshape(find(ones(Nb(whichosc),Nb(whichosc))),[Nb(whichosc),Nb(whichosc)]) );        
                rhovtrace = sum( rhoscs{oo}(:,n_idx) ,2);                         
                for tt=1:Nt
                    % Normalize state
                     rhoscs{oo}(tt,:) = rhoscs{oo}(tt,:) / rhovtrace(tt);
                   %for bb=1:Nb(1)                   
                   %   avgN(tt,oo)  = avgN(tt,oo)  + (bb-1)   * rhoscs{oo}(tt,n_idx(bb));
                   %   avgN2(tt,oo) = avgN2(tt,oo) + (bb-1)^2 * rhoscs{oo}(tt,n_idx(bb));
                   %end
                end

                subplot(lastosc,1,oo)        

                    plot(time,rhoscs{oo}(:,n_idx),'LineWidth',1,'Color',colores(qq)),hold on
                    plot(time, rhovtrace ,'k--', 'LineWidth',1);
                    grid on
                    xlabel('Time(fs)'      ,'interpreter','latex','FontSize',15)
                    ylabel('$<n|\hat{\rho}|n>$','interpreter','latex','FontSize',15)
                    title(sprintf('osc. %i', oo),'interpreter','latex')
                    xlim([tmin tlim])
            end
        end
        
        if 0
        % Mandel parameter
        figure(4);clf;figure(4)

            mandel = zeros(Nt,numel(getreducedosc));
            for tt=1:Nt
            for oo=1:numel(getreducedosc)    
               %mandel(tt,oo) =  (occup2(tt,oo)-occup(tt,oo)^2) / occup(tt,oo) - 1;
                mandel(tt,oo) =  (avgN2(tt,oo)-avgN(tt,oo)^2) / avgN(tt,oo) - 1;
            end
            end

            mandelneg = real(mandel);
            mandelneg(mandelneg>0) = 0;

            coloresMand=['k','r','b','g','y','c','m','k','r','b','g','y','c','m'];
            %coloresMand=[      0    0.4470    0.7410;
            %                    0.8500    0.3250    0.0980;
            %                    0.9290    0.6940    0.1250;
            %                    0.4940    0.1840    0.5560;
            %                    0.4660    0.6740    0.1880;
            %                    0.3010    0.7450    0.9330;
            %                    0.6350    0.0780    0.1840;
            %                    0    0.4470    0.7410;
            %                    0.8500    0.3250    0.0980;
            %                    0.9290    0.6940    0.1250;
            %                    0.4940    0.1840    0.5560;
            %                    0.4660    0.6740    0.1880;
            %                    0.3010    0.7450    0.9330;
            %                    0.6350    0.0780    0.1840];

            for i=1:size(mandelneg,2)     
               %area1= area( time, mandelneg(:,i), 0 ,'FaceColor', coloresMand(i) ,'LineWidth',.002,'LineStyle','-'); hold on
                area1= area( time, mandelneg(:,i), 0 ,'LineWidth',.002,'LineStyle','-'); hold on
               if i==1
                   hold on
               end           
               area1.FaceAlpha=0.3;
            end


            for oo=1:numel(getreducedosc)
                if oo==1 || oo==2
                plot(time,mandel(:,oo),coloresMand(oo),'LineWidth',2); hold on
%                 plot(time,mandel(:,oo),'LineWidth',2); hold on
                else
                plot(time,mandel(:,oo),'LineWidth',0.1,'LineStyle','-','Color',[0.6 0.6 0.6]); hold on    
                end
            end
            
            h1=plot(time,mandel(:,1),'k','LineWidth',1); hold on
            h2=plot(time,mandel(:,2),'r','LineWidth',1);    
            h3=plot(time,mandel(:,3),'b','LineWidth',1);    
            h4=plot(time,mandel(:,4),'f','LineWidth',1);    
            h5=plot(time,mandel(:,5),'y','LineWidth',1);    
            h6=plot(time,mandel(:,5),'y','LineWidth',1);    
            h6=plot(time,mandel(:,5),'y','LineWidth',1);    

            ylabel('$Q(t)$','interpreter','latex','FontSize',15)
            xlabel('t(fs)'      ,'interpreter','latex','FontSize',15)
            title('Mandel~parameter $Q(t)=\frac{<\hat n^2>-<\hat n>^2}{<\hat n>}-1$','interpreter','latex')
            xlim([tmin,tlim])     
            grid on
            
            %darkBackground(fig1,[0.1 0.1 0.1])
            %%darkBackground(fig2,[0.1 0.1 0.1])
            %darkBackground(fig3,[0.1 0.1 0.1])
            %darkBackground(fig4,[0.1 0.1 0.1],[1 0.1 0.1])
        end
            
    end
    
        if qq ~= numel(source)                                  
            disp('Click Me Baby One more Time')
            pause         
        else
        end
    end
    disp('Finished plotting data')
end