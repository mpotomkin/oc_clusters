% IDENTIFIES CLUSTERS AND CLUSTER NUMBERS.
% DATA IS TAKEN FROM FOLDER "numerical_results"

% first clean figures and memory from all previous calculations
clc;
clf;
clear;

Nc = 32; 
M = 101;
number_of_realizations = 20;
dist_char=2.5;
L=20;
d_thetta=pi/12;
thetta=[0:d_thetta:pi];
colors=["red";"blue";"magenta";"green";"yellow";"black";"red";"blue";"magenta";"green";"yellow";"black"];


for realization = 1:number_of_realizations
    ch = sprintf("numerical_results/cells%d.txt",realization);
    trajectories = readmatrix(ch);
   
    Nt=size(trajectories,1)/(Nc);
    
    trajectories_frame_x=trajectories(:,1);
    trajectories_frame_y=trajectories(:,2);
    
    clusters_recorded=0; 
    for t=1:Nt
        x=trajectories(Nc*(t-1)+1:Nc*t,1);
        y=trajectories(Nc*(t-1)+1:Nc*t,2);
    
        v_x=trajectories(Nc*(t-1)+1:Nc*t,3);
        v_y=trajectories(Nc*(t-1)+1:Nc*t,4);
        
        phi=trajectories(Nc*(t-1)+1:Nc*t,5);
 
        is_alive = trajectories(Nc*(t-1)+1:Nc*t,8);
        
    end
    

    % COMPUTING CLUSTER STATISTICS FOR THE LAST FRAME  
    Ncc = sum(is_alive); % the current number of cells 
    % compute distance between cells accounting for periodicity 
    for i=1:Ncc 
        for j=1:Ncc
            r=sqrt((per(x(j),x(i)))^2+(per(y(j),y(i)))^2);
            weight(i,j)=(r<dist_char);
        end
    end
    
    % initializing cell labels, the number of cluster which a given cell
    % belongs to. 
    for i=1:Nc
        label(i)=0;
    end
    
    % assigning actual labels- if weight(i,j)==1 than cells i and j belong
    % to the same cluster
    cluster_counter=0;
    for i=1:Ncc
        if (label(i)==0) 
            cluster_counter=cluster_counter+1;
            label(i)=cluster_counter;
            new=1;
        end
        while (new>0) 
            new=0;
            for j1=1:Ncc
                if (label(j1)==0) 
                    for j2=1:Ncc
                        if (label(j2)==label(i))&&(weight(j1,j2)) 
                            label(j1)=label(i);
                            new=new+1;
                        end
                    end
                end
            end
        end
    end
    
    
    for k=1:Ncc
            cluster(k)=0;
            x_cluster_center(k)=0;
            y_cluster_center(k)=0;
            thetta_orientation(k)=0;
            sum_dist_min(k)=0;
    end
    
    % Compute how many cells in each cluster
    for k=1:Ncc
        for i=1:Ncc
            cluster(k)=cluster(k)+(label(i)==k);
        end
    end

    for k=1:Ncc
        if (cluster(k)>=4)  
            %moving cluster to the center to avoid any problems with periodicity
            nc=0;
            for i=1:Ncc
                if (label(i)==k)
                    nc=nc+1;
                    if (nc==1) 
                        xc(nc)=x(i);
                        yc(nc)=y(i);
                    end
                    if (nc>1) 
                        xc(nc)=L/2+signed_per(x(i),xc(1));
                        yc(nc)=L/2+signed_per(y(i),yc(1));
                    end
                end
            end
            x_cluster_shift(k)=xc(1)-L/2;
            y_cluster_shift(k)=yc(1)-L/2;
            xc(1)=L/2;
            yc(1)=L/2;
        
        
            %computing the center
            for i=1:cluster(k)
                x_cluster_center(k)=x_cluster_center(k)+xc(i);
                y_cluster_center(k)=y_cluster_center(k)+yc(i);
            end
            x_cluster_center(k)=x_cluster_center(k)/cluster(k);
            y_cluster_center(k)=y_cluster_center(k)/cluster(k);
        
            %computing covariance matrices 
            aa=0;
            bb=0;
            cc=0;
            dd=0;
            for i=1:cluster(k)
               j=i;
               aa=aa+(xc(i)-x_cluster_center(k))*(xc(j)-x_cluster_center(k));
               bb=bb+(xc(i)-x_cluster_center(k))*(yc(j)-y_cluster_center(k));
               cc=cc+(yc(i)-y_cluster_center(k))*(xc(j)-x_cluster_center(k));
               dd=dd+(yc(i)-y_cluster_center(k))*(yc(j)-y_cluster_center(k));
            end
            aa=aa/(cluster(k))^2;
            bb=bb/(cluster(k))^2;
            dd=dd/(cluster(k))^2;
                       
            lambda_plus=0.5*(aa+dd+sqrt((aa-dd)^2+4*bb^2));
            lambda_minus=0.5*(aa+dd-sqrt((aa-dd)^2+4*bb^2));
            
            lambda_m=(abs(lambda_plus)>abs(lambda_minus))*(lambda_plus-lambda_minus)+lambda_minus;
            
            thetta0_orientation(k)=atan((lambda_m-aa)/bb);
            
            
            %computing the orientation
            sum_dist=0;
            for i=1:cluster(k)
                sum_dist=sum_dist+abs(cos(thetta(1))*(xc(i)-x_cluster_center(k))+sin(thetta(1))*(yc(i)-y_cluster_center(k)));
            end
            sum_dist_min(k)=sum_dist;
            thetta_orientation(k)=thetta(1);
            for j=2:length(thetta)
                sum_dist=0;
                for i=1:cluster(k)
                    sum_dist=sum_dist+abs(cos(thetta(j))*(xc(i)-x_cluster_center(k))+sin(thetta(j))*(yc(i)-y_cluster_center(k)));
                end
                if (sum_dist<sum_dist_min(k))  
                    sum_dist_min(k)=sum_dist;
                    thetta_orientation(k)=thetta(j);
                end
                
            end
            
            %determining width and length of the cluster
            for i=1:cluster(k)
                wid(i)=abs(cos(thetta_orientation(k))*(xc(i)-x_cluster_center(k))+sin(thetta_orientation(k))*(yc(i)-y_cluster_center(k)));
                len(i)=abs(sin(thetta_orientation(k))*(xc(i)-x_cluster_center(k))-cos(thetta_orientation(k))*(yc(i)-y_cluster_center(k)));
            end
            width_cluster(k)=max(wid(1:cluster(k)));
            length_cluster(k)=max(len(1:cluster(k)));
        
             %determining width and length of the cluster
            for i=1:cluster(k)
                wid(i)=abs(cos(thetta0_orientation(k))*(xc(i)-x_cluster_center(k))+sin(thetta0_orientation(k))*(yc(i)-y_cluster_center(k)));
                len(i)=abs(sin(thetta0_orientation(k))*(xc(i)-x_cluster_center(k))-cos(thetta0_orientation(k))*(yc(i)-y_cluster_center(k)));
            end
            width0_cluster(k)=max(wid(1:cluster(k)));
            length0_cluster(k)=max(len(1:cluster(k)));
            
            x_cluster_center(k)=x_cluster_center(k)+x_cluster_shift(k);
            y_cluster_center(k)=y_cluster_center(k)+y_cluster_shift(k);
            
            %determining neighbors
            for m1=1:Ncc
                for m2=1:Ncc
                    connectivity(m1,m2)=0;
                    if (m1~=m2)
                        distance0(1) = sqrt((x(m1)-x(m2))^2+(y(m1)-y(m2))^2);
                        distance0(2) = sqrt((x(m1)-x(m2)+L)^2+(y(m1)-y(m2))^2);
                        distance0(3) = sqrt((x(m1)-x(m2)-L)^2+(y(m1)-y(m2))^2);
                        distance0(4) = sqrt((x(m1)-x(m2))^2+(y(m1)-y(m2)+L)^2);
                        distance0(5) = sqrt((x(m1)-x(m2))^2+(y(m1)-y(m2)-L)^2);
                        distance0(6) = sqrt((x(m1)-x(m2)+L)^2+(y(m1)-y(m2)+L)^2);
                        distance0(7) = sqrt((x(m1)-x(m2)+L)^2+(y(m1)-y(m2)-L)^2);
                        distance0(8) = sqrt((x(m1)-x(m2)-L)^2+(y(m1)-y(m2)+L)^2);
                        distance0(9) = sqrt((x(m1)-x(m2)-L)^2+(y(m1)-y(m2)-L)^2);
                        distance = min(distance0);
                        if distance < dist_char
                            connectivity(m1,m2)=1; % seems the same as weight(m1,m2)
                        end
                    end
                end
            end 
        end
    end
    
    connectivity_index(1)=0;
    for k=1:length(connectivity_index)
        connectivity_index(k) = 0;
    end

    for k=1:Ncc
        total_connections_in_a_cluster(k)=0;
        connectivity_index(k)=0;
        if (cluster(k)>=4) 
            for m=1:Ncc
                if (label(m)==k) 
                    for m2=1:Ncc
                        total_connections_in_a_cluster(k)=total_connections_in_a_cluster(k)+connectivity(m,m2)/2;
                    end
                end
            end
            connectivity_index(k)=total_connections_in_a_cluster(k)/cluster(k);
        end                  
    end
    %computation of S1 from Zheng paper
    for k=1:Ncc 
        n_m(k)=0;
    end
    for k=1:Ncc
        if (cluster(k)>0)  
            n_m(cluster(k))=n_m(cluster(k))+1;
        end
    end
    
    summ_1=0;
    summ_2=0;
    summ_3=0;%number of non-zero clusters 
    summ_4=0;
    for m=1:Ncc 
        summ_1=summ_1+m^2*n_m(m);
        summ_2=summ_2+m*n_m(m);
        summ_3=summ_3+(cluster(m)>0);
        summ_4=summ_4+cluster(m);
    end    
    S1(realization)=summ_1/(summ_2);
    S2(realization)=summ_4/summ_3;    
    CI(realization)=sum(connectivity_index)/sum(connectivity_index~=0);

    figure(realization)
    clf; 

    for i=1:Ncc 
        if (label(i)<12) 
            color_fun = abs(label(i)/12);
            plot(x(i),y(i),'.','Color',colors(label(i)),'MarkerSize',15); hold on;
        else 
            plot(x(i),y(i),'black.');hold on;
        end
    end

    if (cluster_counter<=100) 
        for k=1:Ncc
            if (cluster(k)>=4) 
                %central
                ll=length0_cluster(k);
                
                x1=x_cluster_center(k)+ll*sin(thetta0_orientation(k));
                x2=x_cluster_center(k)-ll*sin(thetta0_orientation(k));
                y1=y_cluster_center(k)-ll*cos(thetta0_orientation(k));
                y2=y_cluster_center(k)+ll*cos(thetta0_orientation(k));
                plot([x1 x2],[y1 y2],'black','LineWidth',2);
                plot([x1+L x2+L],[y1 y2],'black','LineWidth',2);
                plot([x1-L x2-L],[y1 y2],'black','LineWidth',2);
                plot([x1 x2],[y1+L y2+L],'black','LineWidth',2);
                plot([x1 x2],[y1-L y2-L],'black','LineWidth',2);


                ww=width0_cluster(k);
                x1=x_cluster_center(k)+ww*cos(thetta0_orientation(k));
                x2=x_cluster_center(k)-ww*cos(thetta0_orientation(k));
                y1=y_cluster_center(k)+ww*sin(thetta0_orientation(k));
                y2=y_cluster_center(k)-ww*sin(thetta0_orientation(k));
                plot([x1,x2],[y1,y2],'black','LineWidth',2);
                  plot([x1+L x2+L],[y1 y2],'black','LineWidth',2);
                plot([x1-L x2-L],[y1 y2],'black','LineWidth',2);
                plot([x1 x2],[y1+L y2+L],'black','LineWidth',2);
                plot([x1 x2],[y1-L y2-L],'black','LineWidth',2);

                
                % %central
                % ll=length_cluster(k);
                % 
                % x1=x_cluster_center(k)+ll*sin(thetta_orientation(k));
                % x2=x_cluster_center(k)-ll*sin(thetta_orientation(k));
                % y1=y_cluster_center(k)-ll*cos(thetta_orientation(k));
                % y2=y_cluster_center(k)+ll*cos(thetta_orientation(k));
                % plot([x1 x2],[y1 y2],'red');
                % 
                % 
                % ww=width_cluster(k);
                % x1=x_cluster_center(k)+ww*cos(thetta_orientation(k));
                % x2=x_cluster_center(k)-ww*cos(thetta_orientation(k));
                % y1=y_cluster_center(k)+ww*sin(thetta_orientation(k));
                % y2=y_cluster_center(k)-ww*sin(thetta_orientation(k));
                % plot([x1,x2],[y1,y2],'red');
                % %to the right
                
                
            end
        end
    end
    


    ch = sprintf ("S1= %2.1f, S2= %2.1f, CI=%2.1f%",S1(realization),S2(realization),sum(connectivity_index)/sum(connectivity_index~=0));
    title(ch,'FontSize',15);
    axis([0 L 0 L]);
    grid off
    daspect([1 1 1])
        
       

end

writematrix(S1,'Des1.txt');
writematrix(S2,'Des2.txt');
writematrix(CI,'Deci.txt');
%%
figure(11)
clf
plot([1:length(S1)],S1,'blue','LineWidth',2); hold on;
plot([1:length(S2)],S2,'red','LineWidth',2); hold on;
plot([1:length(CI)],CI,'black','LineWidth',2); hold on;
legend({'S_1','S_2','CI'});
xlabel('Degradation rate');

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
   d(1)=abs(a-b+L);
   d(2)=abs(a-b);
   d(3)=abs(a-b-L);
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
