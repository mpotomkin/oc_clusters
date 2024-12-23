% THIS MATLAB GENERATES VIDEO FROM DATA IN "numerical_results"
% VIDEO SHOWS CELL DYNAMICS ON A DEGRADING SUBSTRATE

% first clean figures and memory from all previous calculations
clc;
clf;
clear;

Nc = 32; 
M = 101;
number_of_realizations = 1;
dist_char=2.5;
L=20;
d_thetta=pi/12;
thetta=[0:d_thetta:pi];
r1 = 14;

for realization = r1:r1%1:number_of_realizations
    ch = sprintf("numerical_results/cells%d.txt",realization);
    trajectories = readmatrix(ch);
    ch = sprintf("numerical_results/substrate_density%d.txt",realization); 
    substrate = readmatrix(ch);

    Nt=size(trajectories,1)/(Nc);
    
    trajectories_frame_x=trajectories(:,1);
    trajectories_frame_y=trajectories(:,2);
    
    clusters_recorded=0; 
    % INITIAL
    for t=1:Nt
        x=trajectories(Nc*(t-1)+1:Nc*t,1);
        y=trajectories(Nc*(t-1)+1:Nc*t,2);
    
        v_x=trajectories(Nc*(t-1)+1:Nc*t,3);
        v_y=trajectories(Nc*(t-1)+1:Nc*t,4);
        
        phi=trajectories(Nc*(t-1)+1:Nc*t,5);
 
        is_alive = trajectories(Nc*(t-1)+1:Nc*t,8);
        
        for i=1:M
            xm(i)=L/M*i;
            for j=1:M        
                substrate_frame_rho(i,j)=substrate((t-1)*M+i,j);
            end
        end

        %%
        %Drawing
        figure(1)
        clf;
        plot([0 0],[0 L],'black');hold on;
        plot([0 L],[L L],'black');hold on;
        plot([L L],[L 0],'black');hold on;
        plot([L 0],[0 0],'black');hold on;
        % subplot(1,2,1);
        % for n=1:Nc
        %     if is_alive(n) == 1
        %         plot(x(n),y(n),'red.','MarkerSize',15); hold on;
        %         plot([x(n) x(n)+0.5*cos(phi(n))],[y(n) y(n)+0.5*sin(phi(n))],'black'); hold on;
        %     end
        % end
        % axis([0 L 0 L]);
        % ch = sprintf ("Time step %d, Cell number %d",t,sum(is_alive));
        % sgtitle(ch);
        % grid off;
        % daspect([1 1 1]);
        % set(gca, 'YDir','reverse')
        %subplot(1,2,2);
        colormap(copper);
        imagesc(xm,xm,substrate_frame_rho'); hold on;
        %colorbar('Ticks',[0,0.5,1],...
        % 'TickLabels',{'\rho=0','\rho=0.5','\rho=1'},'FontSize',15)
        clim([0 1]);
        Z=substrate_frame_rho-3;
        %surf(Z,'FaceAlpha',0.0);
        for n=1:Nc
            if is_alive(n) == 1
                plot(x(n),y(n),'green.','MarkerSize',24); hold on;
                %plot([x(n) max(0,min(L,x(n)+0.5*cos(phi(n))))],[y(n) max(0,min(L,y(n)+0.5*sin(phi(n))))],'green'); hold on;
            end
        end
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'ZTick',[])
        
        axis([0 L 0 L]);
        ch = sprintf ("Time step %d, Cell number %d",t,sum(is_alive));
        %title(ch);
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
        grid off;
        daspect([1 1 1]);
        view(3)
        set(gca,'visible','off')


        saveas(gcf, sprintf("figs/%d.png",t));
    end


end 


function z=per(a,b)
    L=20;
    z=min([abs(a-b),abs(a-b+L),abs(a-b-L)]);
end

function z=corr(x)
    L=20;
    z=min(L,max(0,x));
end

function z=ell_corr(x,y,phi,ell)
    L=20;
    helper=ell;
    if x+ell*cos(phi)>L  
        helper=(L-x)/cos(phi);
    end
    if x+ell*cos(phi)<0 
        helper=-x/cos(phi);
    end
    z=helper;
end

function z=signed_per(a,b)
   L=20;
   d(1)=abs(a-b+L)
   d(2)=abs(a-b)
   d(3)=abs(a-b-L)
   m=d(2);
    indicator=0;
   if d(1)<m 
       indicator=-1;
       m=d(1);
   end
   if d(3)<m 
       indicator=1;
   end
    z=a-b-indicator*L;
end
