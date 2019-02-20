function [gbest,gbestval,fitcount,history,saved]= Diversified_DDSRPSO(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
rand('state',sum(100*clock));
me=Max_Gen;
ps=Particle_Number;
D=Dimension;
cc=[1.49445 1.49445];   %acceleration constants
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=(0.100625/1.5)*(VRmax-VRmin);
%mv1 =(0.100625/1.5)*200;
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;
pos=VRmin+(VRmax-VRmin).*rand(ps,D);
e=feval(fhd,pos',varargin{:});
fitcount=ps;
vel=Vmin+(Vmax-Vmin).*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=e; %initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
dis=abs(gbest-pos);
mean_dis=mean(dis);
gbestrep=repmat(gbest,ps,1);
w_varyfor = floor(0.501*me); 
w_now    = 1.05*ones(ps,1); %Initial Inertia Weight
inertdec = ((1.05-0.5)/w_varyfor)*ones(ps,1); % 

K = 3;
p=1-power(1-1/ps,K);

history = [0, gbestval];
dLow=5e-6;
dHigh=0.25;
pbar=mean(pos);
diag=200;
t1=sqrt(sum((pos-pbar).^2,2));
diver=1/(ps*diag)*sum(t1);
dir=1;
iwt=rand(ps,1);
arch_gbest=[];
arch_gbestval=[];
epsilon=0.05;
saved=0;
used=ps;
old_pos=pos;
II=logical(ones(ps,1));
for i=2:me       
        if(diver<dLow && dir==1)
            arch_gbest=[arch_gbest;gbest];
            arch_gbestval=[arch_gbestval;gbestval];
            dir=-1;
        elseif(diver>dHigh && dir<0)
            dir=1;
        end
        if i>2000
            dir=1;
        end
        
        w_now([1:gbestid-1 gbestid+1:ps]) = w_now([1:gbestid-1 gbestid+1:ps]) - inertdec([1:gbestid-1 gbestid+1:ps]);
        w_now(gbestid) = w_now(gbestid) + inertdec(gbestid);
        
        flg = rand(ps,D)>0.5;
        %flg=dis>repmat(mean_dis,size(dis,1),1);%(rand(ps,D).*
        flg1 = rand(ps,D)>0;
        flg1(gbestid,:) = 0;
        flg(gbestid,:) = 0;
        R1 = rand(ps,D);
        R2 = rand(ps,D);
        
        [~,ii] = sort(pbestval);
         worst = ii(end-2:end);
        flg(worst,:) = 1;
        tvstep = vel(worst,:);
        rr = setdiff(ii,worst);
        bb = setdiff(rr,gbestid);
              
             % determine how many elements is required percent
              numelements = round(0.4*length(bb));
              % get the randomly-selected indices
              indices = randperm(length(bb));
              indices = indices(1:numelements);
              % choose the subset of a you want
              b_sel = bb(indices);

          lPBest(1,:) = median(pbest(rr(2:4),:));
          lPBest(2,:) = median(pbest(rr(2:4),:));
          lPBest(3,:) = median(pbest(rr(2:4),:));
           
        aa=cc(1).*flg1.*R1.*(pbest-pos)+cc(2).*flg.*R2.*(gbestrep-pos);
        vel=repmat(w_now,1,D).*vel+dir*aa; 
 
        if ~isempty(b_sel)
            N = numel(b_sel);
            L=zeros(N,N);
            for s=1:1:N
                L(s,s)=1;
            end
            
            for s = 1:1:N % Each particle (column) informs at most K other at random
                for r=1:1:N
                    if (r~=s)
                        if (alea(0,1)<p)
                            L(s,r) = 1;
                        end
                    end
                end
            end       
        
        for n = 1:1:N  % For each particle (line) ..
            
            %  ...find the best informant g
            MIN = Inf;
            for s=1:1:N
                if (L(s,n) == 1)
                    if pbest(s) < MIN
                        MIN = pbest(s);
                        g_best = s;
                    end
                end
            end
            
            % define a point p' on x-p, beyond p
            p_x_p = pos(n,:) + cc(1)*(pbest(n,:) - pos(n,:));
            % ... define a point g' on x-g, beyond g
            p_x_l = pos(n,:) + cc(2)*(pos(g_best,:) - pos(n,:));
            
            if (g_best == n) % If the best informant is the particle itself, define the gravity center G as the middle of x-p'
                G = 0.5*(pos(n,:) + p_x_p);
            else % Usual  way to define G
                sw = 1/3;
                G = sw*(pos(n,:) + p_x_p + p_x_l);
            end
            
            
            rad = norm(G - pos(n,:)); % radius = Euclidean norm of x-G
            x_p = alea_sphere(D,rad)+ G; % Generate a random point in the hyper-sphere around G (uniform distribution)
            vel(n,:) = w_now(n,:).*vel(n,:) + x_p - pos(n,:); % Update the velocity = w*v(t) + (G-x(t)) + random_vector
            % The result is v(t+1)
         %   x(n,:) = x(n,:) + v(n,:); % Apply the new velocity to the current position. The result is x(t+1)
            
        end
        end
   %     vel(b_sel,:)=repmat(w_now(b_sel),1,D).*vel(b_sel,:)+cc(1).*flg1(b_sel,:).*R1(b_sel,:).*(pbest(b_sel,:)-pos(b_sel,:));
        
          
        vel(worst,:) = repmat(w_now(worst),1,D).*tvstep + cc(1).*R1(worst,:).*flg1(worst,:).*(lPBest-pos(worst,:))+ cc(2)*(R2(worst,:)).*flg(worst,:).*(gbestrep(worst,:)-pos(worst,:));
        vel(rr,:)=max(Vmin(rr,:),min(vel(rr,:),Vmax(rr,:)));
        vel(worst,:)=max(Vmin(worst,:)*1.25,min(vel(worst,:),Vmax(worst,:)*1.25));      
        old_pos(II,:)=pos(II,:);
        pos=pos+vel;
          pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
             +(pos<VRmin).*(VRmin+0.25.*(VRmax-VRmin).*rand(ps,D))+(pos>VRmax).*(VRmax-0.25.*(VRmax-VRmin).*rand(ps,D));
         
        diff=sqrt(sum((old_pos-pos).^2,2));
        II=diff>epsilon;
        used=used+ps;
        saved=saved+(ps-sum(II));
         
         
        e=feval(fhd,pos',varargin{:});
        fitcount=fitcount+ps;      
        tmp=(pbestval<e);
        temp=repmat(tmp',1,D);    
        pbest=temp.*pbest+(1-temp).*pos;
        pbestval=tmp.*pbestval+(1-tmp).*e;%update the pbest
        [gbestval,gbestid]=min(pbestval);

        gbest=pbest(gbestid,:);
        gbestrep=repmat(gbest,ps,1);%update the gbest
        history((size(history,1)+1), :) = [i, gbestval];
        %dis=abs(gbest-pos);
        dis=abs(gbest-pos);
        mean_dis=mean(dis);
        
        pbar=mean(pos);
        t1=sqrt(sum((pos-pbar).^2,2));
        diver=1/((ps)*diag)*sum(t1);
        %ranking mechanism
        [~,~,rank]=unique(pbestval);
        iwt=0.4+0.5*(rank./ps);
end
arch_gbestval=[arch_gbestval;gbestval];
gbestval=min(arch_gbestval);
end
