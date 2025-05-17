

function [tissout]=secomb_oxygen(mesh,bcp,mxx,myy,mzz,alx,aly,alz)

%% geomtry,flow,bcs


allsegments=mesh.LINES;
nnod=length(mesh.POS);
nseg=length(mesh.LINES);
nodname=1:nnod;
segnodname=[allsegments(:,1)';allsegments(:,2)'];
ista=segnodname(1,:);
iend=segnodname(2,:);
segtyp=5*ones(size(segnodname));
diam=1*double(mesh.DIAM);
q=mesh.Q;
qq=abs(q);
hd=mesh.hd;



nnodbc=length(mesh.allbcs); %number of boundary condition nodes

bcnod=mesh.allbcs;
cnode=[(mesh.POS(:,1))';(mesh.POS(:,2))';(mesh.POS(:,3))'];


%% inputs/paramters

slsegdiv=5;
solutefac=1;
flowfac=(1e6)/60;
oxygen=1;
q0fac=1;
% tissparam=[0.000167;10;0];
% tissparam=[1.6667e-04;10;0];
tissparam=[1.666666666666667e-04;10;0];
nsp=1; %number of solutes






lb=125; %outer bound distance


maxl=20; %maximum allowable segment length


v=alx*aly*alz;
vol=v/(mxx*myy*mzz);
req=(vol*0.75/pi)^(1/3);


lam=0.2;
greensverbose=0;


nresis=16;
resisdiam=[4;5;6;8;10;12;14;16;18;20;22;24;28;32;40;50];
resis=[2.45;2;1.74;1.44;1.27;1.16;1.08;1.02;0.974;0.935;0.895;0.853;0.815;0.779;0.75;0.73];
intravascfac=1;
fac=1/4/pi;

lowflowcrit=500;
% hd=0.4*ones(size(hd));
g0=0.2;
g0fac=1;
permsolute=1;

pref=90;
errfac=1e-3; %changed for dark

diffsolute=1;
g0method=2;
nmax=100;
diff1=0.06; %D*alpha
nmaxvessel=2;
nmaxtissue=5;
%%

nodtyp=zeros(1,nnod);

for iseg=1:nseg
    
    

    for inod=1:nnod
            if nodname(inod)==segnodname(1,iseg)
        
            ista(iseg)=inod;
            end
    end
    
    for inod=1:nnod
        
            if nodname(inod)==segnodname(2,iseg)
        
            iend(iseg)=inod;
            end
        
    end
    
    
end


for iseg=1:nseg
    
    inod1=ista(iseg);
    inod2=iend(iseg);
    nodtyp(inod1)=nodtyp(inod1)+1;
    nodtyp(inod2)=nodtyp(inod2)+1;
    
    nodseg(nodtyp(inod1),inod1)=iseg;
    nodseg(nodtyp(inod2),inod2)=iseg;

    nodnod(nodtyp(inod1),inod1)=inod2;
    nodnod(nodtyp(inod2),inod2)=inod1;
    
end

%% start[k][iseg]=coordinates of starting point of segment i

for iseg=1:nseg
    rseg(iseg)=diam(iseg)/2;
    lseg(iseg)=0;
    
    for k=1:3
        
        start(k,iseg)=cnode(k,ista(iseg));
        end1(k,iseg)=cnode(k,iend(iseg));
        ss(k)=end1(k,iseg)-start(k,iseg);
        lseg(iseg)=lseg(iseg)+ss(k)^2;
        
    end
    lseg(iseg)=sqrt(lseg(iseg));
    
    for j=1:3       
        scos(j,iseg)=ss(j)/lseg(iseg);      
    end
    
    for kk=0:15
        t=kk*2*pi/16;
        sintheta=sqrt(1-scos(3,iseg)^2);
        
        if sintheta>0.0001
            cosfi=-scos(2,iseg)/sintheta;
            sinfi= scos(1,iseg)/sintheta;
            
        else
            cosfi=1;
            sinfi=0;
        end
        
        rsta(1,kk+1,iseg)=rseg(iseg)*cosfi*cos(t)-rseg(iseg)*scos(3,iseg)*sinfi*sin(t)+start(1,iseg);
        rsta(2,kk+1,iseg)=rseg(iseg)*sinfi*cos(t)+rseg(iseg)*scos(3,iseg)*cosfi*sin(t)+start(2,iseg);
        rsta(3,kk+1,iseg)=rseg(iseg)*sintheta*sin(t)+start(3,iseg);
        
        
        rend(1,kk+1,iseg)=rseg(iseg)*cosfi*cos(t)-rseg(iseg)*scos(3,iseg)*sinfi*sin(t)+end1(1,iseg);
        rend(2,kk+1,iseg)=rseg(iseg)*sinfi*cos(t)+rseg(iseg)*scos(3,iseg)*cosfi*sin(t)+end1(2,iseg);
        rend(3,kk+1,iseg)=rseg(iseg)*sintheta*sin(t)+end1(3,iseg);
        
        
    end
    
    
    
end


%% subdivide segments into smaller segments

nnv=0;

for iseg=1:nseg
    
    m=floor(lseg(iseg)/maxl+1);
    nspoint(iseg)=m;
    istart(iseg)=nnv+1;
    ds(iseg)=lseg(iseg)/m;
    nnv=nnv+m;
    
    
end



%% compute coordinates of tissue points

delx=alx/mxx;
dely=aly/myy;
delz=alz/mzz;

for i=1:mxx   
    axt(i)=(i-0.5)*delx;   
end

for i=1:myy   
    ayt(i)=(i-0.5)*dely;   
end

for i=1:mzz   
    azt(i)=(i-0.5)*delz;   
end


%% outboun

for i=1:mxx
    for j=1:myy
        for k=1:mzz
            nbou(i,j,k)=1;
        end        
    end    
end

am=6;

for a=0:am
    for b=-am:am
        for c=-am:am
            
            if a~=0 || b~=0 || c~=0
                aa=0;
                bb=0;
                cc=0;
                
                if a~=0
                    aa=1/a;
                end
                
                if b~=0
                    bb=1/b;
                end
                
                if c~=0
                    cc=1/c;
                end
                abc=sqrt((aa^2)+(bb^2)+(cc^2));
                ddmin=1e8;
                ddmax=-1e8;
                
                
                for iseg=1:nseg
                    
                    dd=(aa*start(1,iseg)+bb*start(2,iseg)+cc*start(3,iseg))/abc;
                    ddmin=min([dd-rseg(iseg)-lb,ddmin]);
                    ddmax=max([dd+rseg(iseg)+lb,ddmax]);
                    
                    dd=(aa*end1(1,iseg)+bb*end1(2,iseg)+cc*end1(3,iseg))/abc;
                    ddmin=min([dd-rseg(iseg)-lb,ddmin]);
                    ddmax=max([dd+rseg(iseg)+lb,ddmax]);
                end
                
                for i=1:mxx
                    for j=1:myy 
                        for k=1:mzz
                            t=(aa*axt(i)+bb*ayt(j)+cc*azt(k))/abc;
                            
                            if t>ddmax||t<ddmin
                                nbou(i,j,k)=0;
                            end
                            
                            
                            
                        end
                    end

                end
                
                
                
            end
            
        end
        
    end
    
end


nnt=0;
for i=1:mxx
    for j=1:myy
        for k=1:mzz
            
            if nbou(i,j,k)==1
          nnt=nnt+1;
          nbou(i,j,k)=nnt;
            end
        end
        
    end
    
end


%%
%%%%%%% greens %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% set up main seg
%mainseg each point on all subsegments with original segment label
mainseg=zeros(1,nnv);

for iseg=1:nseg
    
    for i=0:nspoint(iseg)-1
        
        mainseg(istart(iseg)+i)=iseg;
    end
    
end

%% identify vessel coordinates (ax) on subsegments

% ax ;  col=subsegment point, row=x,y,z

for i=1:nnv
    
    iseg=mainseg(i);
    isn=i-istart(iseg);
    
    for j=1:3
        ax(j,i)=start(j,iseg)+scos(j,iseg)*ds(iseg)*(isn+0.5);
        
    end
    
end

%% index tissue points


bello2=zeros(size(nbou)); %1992


for i=1:mxx
    for j=1:myy
        
        for k=1:mzz
            nt=nbou(i,j,k);
            pr=bello2(i,j,k);
            if nt>0
                
                tisspoints(1,nt)=i;
                tisspoints(2,nt)=j;
                tisspoints(3,nt)=k;

                prpoints(1,nt)=0;
                prpoints(2,nt)=0;
                prpoints(3,nt)=0;
                            if pr>0 
                                 prpoints(1,nt)=1;
                                 prpoints(2,nt)=1;
                                 prpoints(3,nt)=1;
                            end

            end



            



            
            
        end
    end
    
    
end


%% calculate distance of tp to nearest vessel

dtave=0;

for itp=1:nnt
    i=tisspoints(1,itp);
    j=tisspoints(2,itp);
    k=tisspoints(3,itp);
    dtmin(itp)=1e6;
    
    for jseg=1:nseg
        
        x11=(axt(i)-start(1,jseg))*scos(2,jseg)-(ayt(j)-start(2,jseg))*scos(1,jseg);
        x22=(ayt(j)-start(2,jseg))*scos(3,jseg)-(azt(k)-start(3,jseg))*scos(2,jseg);
        x33=(azt(k)-start(3,jseg))*scos(1,jseg)-(axt(i)-start(1,jseg))*scos(3,jseg);
        
        
        
        
        disp2=(x11^2)+(x22^2)+(x33^2);
        ds2=((axt(i)-start(1,jseg))^2)+((ayt(j)-start(2,jseg))^2)+((azt(k)-start(3,jseg))^2);
        de2=((axt(i)-end1(1,jseg))^2)+((ayt(j)-end1(2,jseg))^2)+((azt(k)-end1(3,jseg))^2);
        
        
        hi(itp,jseg)=x33;
            
        
        if (max([ds2,de2])-disp2)>(lseg(jseg)^2)
            
           d=max([sqrt(min([ds2,de2]))-rseg(jseg),0]);

        else
            d=max([sqrt(disp2)-rseg(jseg),0]);

        end
        
        
        
        if d<dtmin(itp)
            dtmin(itp)=d;
        end
        
    end
    dtave=dtave+dtmin(itp);
end
%%
dtave=dtave/nnt;
vdom=nnt*vol;
tlength=0;
tlengthq=0;
tlengthqhd=0;


qdata=q;
for iseg=1:nseg
    q(iseg)=qdata(iseg)*q0fac;
    qq(iseg)=abs(q(iseg));
    tlength=tlength+lseg(iseg);
    tlengthq=tlengthq+lseg(iseg)*qq(iseg);
    tlengthqhd=tlengthqhd+lseg(iseg)*qq(iseg)*hd(iseg);
end

den=sqrt(vdom/tlength);

for isp=1:nsp
    for iseg=1:nseg
        gamma1(iseg,isp)=0;
        if nresis(isp) ~=0
            gamma1(iseg,isp)=resis(1,isp);
            
            for j=2:nresis(isp)
                if diam(iseg)<=resisdiam(j,isp)&&diam(iseg)>resisdiam(j-1,isp)
                    gamma1(iseg,isp)=resis(j-1,isp)+(resis(j,isp)-resis(j-1,isp))*(diam(iseg)-resisdiam(j-1,isp))/(resisdiam(j,isp)-resisdiam(j-1,isp));
                end       
            end
            
           if diam(iseg)>resisdiam(nresis(isp),isp)
               gamma1(iseg,isp)=resis(nresis(isp),isp);
           end
           if oxygen(isp)~=1
               gamma1(iseg,isp)=gamma1(iseg,isp)/pi/diam(iseg);
           end
            
           gamma1(iseg,isp)=gamma1(iseg,isp)*intravascfac(isp);
            
        end

    end    
end


%% gvv


for i=1:nnv
    iseg=mainseg(i);
    
    for j=1:nnv
        jseg=mainseg(j);
        dist=sqrt(((ax(1,j)-ax(1,i))^2)+((ax(2,j)-ax(2,i))^2)+((ax(3,j)-ax(3,i))^2));
        
        if dist<max([sqrt(ds(iseg)*rseg(iseg)),sqrt(ds(jseg)*rseg(jseg))])
            dsmax=max([ds(iseg),ds(jseg)]);
            rsegmax=max([rseg(iseg),rseg(jseg)]);
            gvarfac=0.6*exp(-0.45*dsmax/rsegmax);
            
            if iseg~=jseg
                dist=rsegmax;   
            end
            gvv(i,j)=(1.298/(1+0.297*(dsmax/rsegmax)^0.838)-gvarfac*((dist/rsegmax)^2))*fac/rsegmax;
        else
            
            gvv(i,j)=fac/dist;
        end
        
    end
    
    
end


%% tissue vessel vessel tissue


gtt=1.2*fac/req;

for jx=1:mxx
    for jy=1:myy
        
        for jz=1:mzz
            
            dist=sqrt(((axt(1)-axt(jx))^2)+((ayt(1)-ayt(jy))^2)+((azt(1)-azt(jz))^2));
            
            if jx*jy*jz~=1
                
                dtt(jx,jy,jz)=fac/dist;
                
            else
                
                dtt(jx,jy,jz)=gtt;
                
            end
            
        end
        
    end
    
end


qhdcrit=0;
for isp=1:nsp
    
    qhdcrit=lowflowcrit*tissparam(1,isp);
    
end

%detect low flow*hd
for iseg=1:nseg
    lowflow(iseg)=0;
    if qq(iseg)*(hd(iseg)+0.01)<qhdcrit
        lowflow(iseg)=1;
    end
end


%% initgreens

for isp=1:nsp
    pinit(isp)=g0(isp);
end


[mtiss,mptiss]=tissrate3(tissparam,nsp,pinit,0); %eq(2)

for isp=1:nsp
    mptissref(isp)=mptiss(isp);
    qtsum(isp)=0;
    for itp=1:nnt
        qt(itp,isp)=mtiss(isp)*vol;
        pt(itp,isp)=pinit(isp);
        qtsum(isp)=qtsum(isp)+qt(itp,isp);
    end
    
    for i=1:nnv
        qv(i,isp)=0;
        pv(i,isp)=0;
        
        if permsolute(isp)==1
            if oxygen(isp)==1
                qv(i,isp)=-qtsum(isp)*ds(mainseg(i))*qq(mainseg(i))*hd(mainseg(i))/tlengthqhd;
                if lowflow(mainseg(i)==1)
                    qv(i,isp)=0;
                end
            else
                qv(i,isp)=-qtsum(isp)*ds(mainseg(i))*qq(mainseg(i))/tlengthq;
            end
            
            pv(i,isp)=pinit(isp);
        end
        
        
        
        
    end
    
end

for isp=1:nsp
    pinit(isp)=0;
end

for isp=1:nsp
    pinit(isp)=pref(isp);
end

for isp=1:nsp
    pinit(isp)=pref(isp);
    
    [mtiss,mptiss]=tissrate3(tissparam,nsp,pinit,0);
    
    epstissue(isp)=max([abs(mtiss(isp))*vol*errfac,0.001]);
    epsvessel(isp)=nnt*epstissue(isp)/nnv;
    eps(isp)=pref(isp)*errfac;
    pinit(isp)=0;

    
end



%% putrank

            for inod=1:nnod
                nodtyp(inod)=0;
                nodout(inod)=0;
            end
            nsegfl=0;
            for iseg=1:nseg
                if q(iseg)>=0
                    nod1=ista(iseg);
                    nod2=iend(iseg);
                else
                    nod1=iend(iseg);
                    nod2=ista(iseg);
                end
                
                nodtyp(nod1)=nodtyp(nod1)+1;
                nodseg(nodtyp(nod1),nod1)=iseg;
                nodnod(nodtyp(nod1),nod1)=nod2;
                nodout(nod1)=nodout(nod1)+1;
                nsegfl=nsegfl+1;
                
                
            end
            
            
            for iseg=1:nseg
                if q(iseg)>=0
                    nod1=ista(iseg);
                    nod2=iend(iseg);
                else
                    nod1=iend(iseg);
                    nod2=ista(iseg);
                end
                nodtyp(nod2)=nodtyp(nod2)+1;
                nodseg(nodtyp(nod2),nod2)=iseg;
                nodnod(nodtyp(nod2),nod2)=nod1;
                
            end
            
            nnodfl=0;
            for inod=1:nnod
                nk(inod)=0;
                if nodtyp(inod)==1 && nodout(inod)==1
                    nnodfl=nnodfl+1;
                    nk(inod)=1;
                    nodrank(nnodfl)=inod;
                end
                
            end
            
            
            
            
            
            
            
            
            flag=1;
            while flag==1
                flag=0;
            for  inod=1:nnod
                
                if nk(inod)==0 && nodtyp(inod)>0
                    
                    
                    for j=nodout(inod)+1:nodtyp(inod)
                        jseg=nodseg(j,inod);
                        if (inod==iend(jseg) && (nk(ista(jseg))==0 || q(jseg)<=0))
                            sl=1;
                            break
                        else
                                                    if (inod==ista(jseg) && (nk(iend(jseg))==0 || q(jseg)>=0))
                                                        sl=1;
                                                        break

                                                    end
                            
                            
                            
                        end

                    end
                    if sl==1
                        sl=0;
                    else
                    nnodfl=nnodfl+1;
                    nk(inod)=1;
                    nodrank(nnodfl)=inod;  
                    flag=1;
                    end
                

                end
                
            end
            end




%%
w_oval=zeros(nnv,nnt);
% % figure
% % hold on
% for i=1:nnv
% 
%     disp(append('calc oval weight ',num2str(i),' out of ',num2str(nnv)))
%     for itp=1:nnt
%         if ax(3,i)<min(ax(3,:))+1
%             
%         n=100;
%         rcirc=0.5*diam(mainseg(i));
%         rad_oval=(rcirc^2)/(thickness_cc/2);
%         if i<nnv
%         curve=[[ax(1,i),ax(1,i+1)];[ax(2,i),ax(2,i+1)];[ax(3,i),ax(3,i+1)]]; %ax(1,i),ax(1,i+1)
%         centre=[ax(1,i),ax(2,i),ax(3,i)];%ax(1,i)
%         else
%         curve=[[ax(1,i),ax(1,i-1)];[ax(2,i),ax(2,i-1)];[ax(3,i),ax(3,i-1)]]; %ax(1,i),ax(1,i+1)
%         centre=[ax(1,i),ax(2,i),ax(3,i)];%ax(1,i)    
%         end
%         [x_oval,y_oval,z_oval]=tubeplot2(curve,[(thickness_cc/2),rad_oval],n);%a =diam/2 %b=thickness/2
%         xxx_oval=x_oval(1:n,2);
%         yyy_oval=y_oval(1:n,2);
%         zzz_oval=z_oval(1:n,2);
% 
%         if max(zzz_oval)>min(ax(3,:))+thickness_cc/2
% 
%         [x_oval,y_oval,z_oval]=tubeplot2(curve,[rad_oval,(thickness_cc/2)],n);%a =diam/2 %b=thickness/2
%         xxx_oval=x_oval(1:n,2);
%         yyy_oval=y_oval(1:n,2);
%         zzz_oval=z_oval(1:n,2);
% 
%         end
% 
% 
%         w=(sqrt(((axt(tisspoints(1,itp))-xxx_oval).^2)+((ayt(tisspoints(2,itp))-yyy_oval).^2)+((azt(tisspoints(3,itp))-zzz_oval).^2)));
%         [~,closest_t_idx]=min(sqrt(((axt(tisspoints(1,itp))-xxx_oval).^2)+((ayt(tisspoints(2,itp))-yyy_oval).^2)+((azt(tisspoints(3,itp))-zzz_oval).^2)));
%         thispoint=w(closest_t_idx);
% %         w_oval(i,itp)=thispoint;
%                 bruh=sort(w);
%         w_oval(i,itp)=mean(bruh(1:86));
% 
% %         [x_circ,y_circ,z_circ]=tubeplot(curve,rcirc,n);
% %         xxx_circ=x_circ(1:n,2);
% %         yyy_circ=y_circ(1:n,2);
% %         zzz_circ=z_circ(1:n,2);
% %         w=(sqrt(((axt(tisspoints(1,itp))-xxx_oval).^2)+((ayt(tisspoints(2,itp))-yyy_oval).^2)+((azt(tisspoints(3,itp))-zzz_oval).^2)))-...
% %         (sqrt(((axt(tisspoints(1,itp))-xxx_circ).^2)+((ayt(tisspoints(2,itp))-yyy_circ).^2)+((azt(tisspoints(3,itp))-zzz_circ).^2)));
% %         [~,closest_t_idx]=min(sqrt(((axt(tisspoints(1,itp))-xxx_oval).^2)+((ayt(tisspoints(2,itp))-yyy_oval).^2)+((azt(tisspoints(3,itp))-zzz_oval).^2)));
% %         thispoint=w(closest_t_idx);
% %         w_oval(i,itp)=0.1*thispoint;
% 
% 
% 
% 
% 
% 
% 
% 
%      
%         else
% 
%         n=100;
%         rcirc=0.5*diam(mainseg(i));
% 
%                 if i<nnv
%                 curve=[[ax(1,i),ax(1,i+1)];[ax(2,i),ax(2,i+1)];[ax(3,i),ax(3,i+1)]]; %ax(1,i),ax(1,i+1)
%                 centre=[ax(1,i),ax(2,i),ax(3,i)];%ax(1,i)
%                 else
%                 curve=[[ax(1,i),ax(1,i-1)];[ax(2,i),ax(2,i-1)];[ax(3,i),ax(3,i-1)]]; %ax(1,i),ax(1,i+1)
%                 centre=[ax(1,i),ax(2,i),ax(3,i)];%ax(1,i)    
%                 end
% 
%         [x_circ,y_circ,z_circ]=tubeplot(curve,rcirc,n);
% 
% 
% 
%         xxx_circ=x_circ(1:n,2);
%         yyy_circ=y_circ(1:n,2);
%         zzz_circ=z_circ(1:n,2);
% 
%         w=(sqrt(((axt(tisspoints(1,itp))-xxx_circ).^2)+((ayt(tisspoints(2,itp))-yyy_circ).^2)+((azt(tisspoints(3,itp))-zzz_circ).^2)));
% %         [~,closest_t_idx]=min(sqrt(((axt(tisspoints(1,itp))-xxx_circ).^2)+((ayt(tisspoints(2,itp))-yyy_circ).^2)+((azt(tisspoints(3,itp))-zzz_circ).^2)));
% %         thispoint=w(closest_t_idx);
% %         w_oval(i,itp)=thispoint;
%                 bruh=sort(w);
%         w_oval(i,itp)=mean(bruh(1:91));
% %         w_oval(i,itp)=0;
%         end
% 
%     end
% %        plot3(xxx_oval,yyy_oval,zzz_oval,'r-')
% % xlabel('x')
% % ylabel('y')
% % zlabel('z')
% % axis equal
% 
% end





%% start of main loop

for kmain=1:200
    disp(append('kmain=',num2str(kmain)))
    
    for isp=1:nsp
        for itp=1:nnt
            ptprev(itp,isp)=pt(itp,isp);
        end
        if permsolute(isp)==1
            for i=1:nnv
                pvprev(i,isp)=pv(i,isp);
                
            end
        end
        
        g0old(isp)=g0(isp);
    end
 
    %start of vessel loop%
    %compute contribution pvt from tissue source strength qt%
    
    for i=1:nnv
        for isp=1:nsp
            if diffsolute(isp)==1
                if g0method==1 && permsolute(isp)==1
                    pvt(i,isp)=0;
                else
                    pvt(i,isp)=g0(isp);
                end
            end
        end
        
        for itp=1:nnt %HERE

            %iif ax(3,i) is less than 30
            %dist = dist+(oval-circle)
            %else
        
            dist=sqrt(((ax(1,i)-axt(tisspoints(1,itp)))^2)+((ax(2,i)-ayt(tisspoints(2,itp)))^2)+((ax(3,i)-azt(tisspoints(3,itp)))^2));
%              dist=w_oval(i,itp);
            if dist<=req
                gvt=fac*(1.5-0.5*((dist/req)^2))/req;
            else
                gvt=fac/(dist);
            end
            
           for isp=1:nsp
               if permsolute(isp)==1
                   pvt(i,isp)=pvt(i,isp)+gvt/diff1(isp)*qt(itp,isp);
               end
           end
            
            
        end     
    end
    
    %compute blood solute levels and p02

    for kvessel=1:nmaxvessel
        convflagv=1;
        for isp=1:nsp
            qvsum(isp)=0;
            dqvsumdg0(isp)=0;
            
            if permsolute(isp)==1
                ineg=0;
                ihigh=0;
                %convect
                 [al,cv,segc,isegk,isegkk]=convect(isp,nseg,nnv,nnodbc,bcnod,nodout,oxygen,bcp,solutefac,flowfac,nnodfl,nodrank,nodtyp,nodseg,qq,hd,istart,nspoint,qv,q);       
                 
                 for i=1:nnv
                     iseg=mainseg(i);
                     qvprev(i,isp)=qv(i,isp);
                     if oxygen(isp)==1
                        if lowflow(iseg)~=1
                         if cv(i,isp)<0
                             ineg=ineg+1;
                         end
                         if cv(i,isp)>bloodconc(150,hd(iseg))
                             ihigh=ihigh+1;
                         end
                         [pv(i,isp),dcdp(i,isp)]=blood(cv(i,isp),hd(iseg));

                        end
                     else
                         pv(i,isp)=cv(i,isp);
                         dcdp(i,isp)=1;
                     end  
                 end
                 
                 %create linear system for solving
                 for i=1:nnv
                     iseg=mainseg(i);
                     rhs(i)=pv(i,isp)-pvt(i,isp);
                     for j=1:nnv
                         jseg=mainseg(j);
                         mat(i,j)=gvv(i,j)/diff1(isp)+al(i,j)/dcdp(i,isp)/qq(iseg)/flowfac;
                         if i==j
                             mat(i,j)=mat(i,j)+gamma1(iseg,isp)/ds(iseg);
                         end
                         rhs(i)=rhs(i)+al(i,j)*qv(j,isp)/dcdp(i,isp)/qq(iseg)/flowfac;
                         if oxygen(isp)==1 && lowflow(mainseg(i))==1
                             if i==j
                                 mat(i,j)=1;
                             else
                                 mat(i,j)=0;
                             end
                             
                         end
                         
                     end
                     if oxygen(isp)==1 &&lowflow(iseg)==1
                         rhs(i)=qvprev(i,isp);
                     end
                 end
                 
                 
                 for i=1:nnv
                     rhsl(i)=rhs(i);
                     matx(i)=qv(i,isp);
                 end
                 
                 
                 matx=bicgstab(mat,rhsl,matx,nnv,1e-4,2000);
                 
                 for i=1:nnv
                     qv(i,isp)=matx(i);
                     qvsum(isp)=qvsum(isp)+qv(i,isp);
                 end
                 for i=1:nnv
                     rhsl(i)=-1;
                     matx(i)=0;
                 end
                 
                 matx=bicgstab(mat,rhsl,matx,nnv,1e-4,2000);
                 
                 for i=1:nnv
                     if oxygen(isp)~=1 || lowflow(mainseg(i))~=1
                         dqvsumdg0(isp)=dqvsumdg0(isp)+matx(i);
                     end
                 end
                 
                 
                 for i=1:nnv
                     iseg=mainseg(i);
                     if oxygen(isp)==1 && lowflow(iseg)==1
                         for j=1:3
                             x(j)=ax(j,i)-0.5*scos(j,iseg)*ds(iseg);
                         end
                         %eval
                         p=eval123(slsegdiv,req,x,nnt,nnv,nsp,mainseg,tisspoints,permsolute,diffsolute,fac,axt,ayt,azt,ds,diff1,g0,qt,scos,qv,ax);
                         p(isp)=max([p(isp),0]);
                         pv(i,isp)=p(isp)/2;
                         qvtemp(i)=q(iseg)*flowfac*bloodconc(p(isp),hd(iseg));
                         for j=1:3
                             x(j)=ax(j,i)+0.5*scos(j,iseg)*ds(iseg);
                         end
                         p=eval123(slsegdiv,req,x,nnt,nnv,nsp,mainseg,tisspoints,permsolute,diffsolute,fac,axt,ayt,azt,ds,diff1,g0,qt,scos,qv,ax);
                         pv(i,isp)=pv(i,isp)+p(isp)/2;
                         qvtemp(i)=qvtemp(i)-q(iseg)*flowfac*bloodconc(p(isp),hd(iseg));
                         cv(i,isp)=bloodconc(pv(i,isp),hd(iseg));
                     end
                 end
                 for i=1:nnv
                     if oxygen(isp)==1 && lowflow(mainseg(i))==1
                         qv(i,isp)=0.5*qvtemp(i)+0.5*qvprev(i,isp);                        
                     end            
                 end
                 errvessel(isp)=0;
                 imaxerr=0;
                 errvesselcount(isp)=0;
                 for i=1:nnv
                     dif=qv(i,isp)-qvprev(i,isp);
                     if qv(i,isp)~=0
                         dif=dif*min([1,epsvessel(isp)/errfac/abs(qv(i,isp))]);
                     end
                     if abs(dif)>=errvessel(isp)
                         imaxerrvessel(isp)=mainseg(i);
                         errvessel(isp)=abs(dif);
                     end
                     if abs(dif)>epsvessel(isp)
                         errvesselcount(isp)=errvesselcount(isp)+1;
                     end
                 end
                 if errvesselcount(isp)>0
                     convflagv=0;
                 end    
            end

        end
        %here
        if convflagv==1
            disp('vessels converged')  
            break
        else
            disp(append(num2str(errvesselcount),'vessels not converged','-error-',num2str(errvessel)))
        
            
            
            continue
        end
    end
    
    %end of vessel loop
    %start of tissue loop
    for itp=1:nnt
        for isp=1:nsp
            ptv(itp,isp)=0;
        end
        for i=1:nnv
            dist=sqrt(((ax(1,i)-axt(tisspoints(1,itp)))^2)+((ax(2,i)-ayt(tisspoints(2,itp)))^2)+((ax(3,i)-azt(tisspoints(3,itp)))^2));
%             dist=w_oval(i,itp);
            if dist<=req
                gtv=fac*(1.5-0.5*((dist/req)^2))/req;
            else
                gtv=fac/(dist);
            end
            for isp=1:nsp
                if permsolute(isp)==1
                    ptv(itp,isp)=ptv(itp,isp)+gtv/diff1(isp)*qv(i,isp);
                end
            end
        end
    end
    for isp=1:nsp
        qvfac(isp)=1;
    end
    
    for ktissue=1:nmaxtissue
        convflagt=1;
        for isp=1:nsp
            qtsum(isp)=0;
            errtissue(isp)=0;
            dqtsumdg0(isp)=0;
            errtissuecount(isp)=0;
        end
        
        for itp=1:nnt
            ix=tisspoints(1,itp);
            iy=tisspoints(2,itp);
            iz=tisspoints(3,itp);
        
        for isp=1:nsp
            ptt(isp)=0;
        end
        for jtp=1:nnt
            jx=tisspoints(1,jtp);
            jy=tisspoints(2,jtp);
            jz=tisspoints(3,jtp);
            ixdiff=abs(ix-jx)+1;
            iydiff=abs(iy-jy)+1;
            izdiff=abs(iz-jz)+1;
            for isp=1:nsp
                if diffsolute(isp)==1
                    ptt(isp)=ptt(isp)+dtt(ixdiff,iydiff,izdiff)*qt(jtp,isp);
                end
            end
        end
        
        for isp=1:nsp
            ptpt(isp)=pt(itp,isp);
        end
        for isp=1:nsp
            if diffsolute(isp)==1
                pt(itp,isp)=(1-lam)*pt(itp,isp)+lam*(ptv(itp,isp)*qvfac(isp)+g0(isp)+ptt(isp)/diff1(isp));
            
            else
                [mtiss,mptiss]=tissrate3(tissparam,nsp,ptpt,prpoints(1,itp));
                pt(itp,isp)=pt(itp,isp)-mtiss(isp)/mptiss(isp);
            end
            ptpt(isp)=pt(itp,isp);
        end
        [mtiss,mptiss]=tissrate3(tissparam,nsp,ptpt,prpoints(1,itp));
        for isp=1:nsp
            dif=mtiss(isp)*vol-qt(itp,isp);
            qt(itp,isp)=mtiss(isp)*vol;
            qtsum(isp)=qtsum(isp)+qt(itp,isp);
            if diffsolute(isp)==1
                dqtsumdg0(isp)=dqtsumdg0(isp)+mptiss(isp)*vol;
            end
            if abs(dif)>errtissue(isp)
                errtissue(isp)=abs(dif);
                imaxerrtissue(isp)=itp;
            end
            if abs(dif)>epstissue(isp)
                errtissuecount(isp)=errtissuecount(isp)+1;
            end
        end
        end
        %632
        
        for isp=1:nsp
            if errtissuecount(isp)>0
                convflagt=0;
            end
        end
        
        
        if kmain>1 && convflagt==1
            break
        end


    end
    
    kvessel=min([kvessel,nmaxvessel]);
    ktissue=min([ktissue,nmaxtissue]);
    %end of tissue loop
    for isp=1:nsp
        g0facnew(isp)=0;
    end
    
    for itp=1:nnt
        for isp=1:nsp
            ptpt(isp)=pt(itp,isp);
        end
        [mtiss,mptiss]=tissrate3(tissparam,nsp,ptpt,prpoints(1,itp));
        
        ix=tisspoints(1,itp);
        iy=tisspoints(2,itp);
        iz=tisspoints(3,itp);
        
        for jtp=1:nnt
            jx=tisspoints(1,jtp);
            jy=tisspoints(2,jtp);
            jz=tisspoints(3,jtp);
            ixdiff=abs(ix-jx)+1;
            iydiff=abs(iy-jy)+1;
            izdiff=abs(iz-jz)+1;
            
            for isp=1:nsp
                if diffsolute(isp)==1
                    g0facnew(isp)=g0facnew(isp)+dtt(ixdiff,iydiff,izdiff)/diff1(isp)*mptiss(isp)*vol;
                    
                end
            end
            
        end
    end
    
    for isp=1:nsp
        if diffsolute(isp)==1 && (g0method==2 || permsolute(isp)==0)
            g0facnew(isp)=1/(1-g0facnew(isp)/nnt);
            dqsumdg0=min([dqvsumdg0,0])+min([dqtsumdg0(isp),0])*g0facnew(isp);
            if abs(dqsumdg0)>1e-6
                dif=(qvsum(isp)+qtsum(isp))/dqsumdg0*g0fac(isp);
                g0(isp)=g0(isp)-dif;
            end
        end
    end
    
    
    %test for convergance
    convflag=1;
    for isp=1:nsp
        err=0;
        imaxerr=0;
        if permsolute(isp)==1
            for i=1:nnv
                dif=abs(pv(i,isp)-pvprev(i,isp))/eps(isp);
                difstore{i,kmain}=dif;
                
                if dif>err
                    imaxerr=mainseg(i);
                    err=dif;
                end
            end
        end
        errvessel(isp)=err;
        imaxerrvessel(isp)=imaxerr;
        
        err=0;
        imaxerr=0;
        for itp=1:nnt
            dif=abs(pt(itp,isp)-ptprev(itp,isp))/eps(isp);
            if dif>err
                imaxerr=itp;
                err=dif;
            end
        end
        errtissue(isp)=err;
        imaxerrtissue(isp)=imaxerr;
        if errvessel(isp)>err
            imaxerr=imaxerrvessel(isp);
            err=errvessel(isp);
        else
            imaxerr=-imaxerr;
        end
        dif=abs(g0(isp)-g0old(isp))/eps(isp);
        if dif>err
            imaxerr=0;
            err=dif;
        end
        
        if err>1
            convflag=0;
        end 

        disp(append('conv err=',num2str(err)))
    end
    if convflag && convflagv && convflagt ==1
        
        break
    end
   
end


for isp=1:nsp
    if diffsolute(isp)==1
        if permsolute(isp)==1
        end
        
        for i=1:nnv
            pev(i,isp)=pvt(i,isp);
            if permsolute(isp)==1
                for j=1:nnv
                    pev(i,isp)=pev(i,isp)+gvv(i,j)*qv(j,isp)/diff1(isp);
                end
                
            end
        end
        for iseg=1:nseg
            if segtyp(iseg)==4 || segtyp(iseg)==5
                qvseg(iseg,isp)=0;
                pevseg(iseg,isp)=0;
                pvseg(iseg,isp)=0;
            end
        end
        
        for i=1:nnv
            iseg=mainseg(i);
            if permsolute(isp)==1
                qvseg(iseg,isp)=qvseg(iseg,isp)+qv(i,isp);
                pevseg(iseg,isp)=pevseg(iseg,isp)+pev(i,isp)/nspoint(iseg);
                pvseg(iseg,isp)=pvseg(iseg,isp)+pv(i,isp)/nspoint(iseg);
            end
        end
        
    end 
end

for isp=1:nsp
    if permsolute(isp)==1
        pmax(isp)=-1e-8;
        pmeanv(isp)=0;
        pmin(isp)=1e8;
        for i=1:nnv
            pmeanv(isp)=pmeanv(isp)+pv(i,isp);
            pmax(isp)=max([pv(i,isp),pmax(isp)]);
            pmin(isp)=min([pv(i,isp),pmin(isp)]);
        end
        pmeanv(isp)=pmeanv(isp)/nnv;
        
    end
end


for isp=1:nsp
    if permsolute(isp)==1
        pmax(isp)=-1e-8;
        pmeant(isp)=0;
        pmin(isp)=1e8;
        for itp=1:nnt
            pmeant(isp)=pmeant(isp)+pt(itp,isp);
            pmax(isp)=max([pt(itp,isp),pmax(isp)]);
            pmin(isp)=min([pt(itp,isp),pmin(isp)]);
        end
        pmeant(isp)=pmeant(isp)/nnt;
        
    end
end


%%
hi=pt;
for iii=1:length(axt)
    for jjj=1:length(ayt) 
        for kkk=1:length(azt)
            if nbou(iii,jjj,kkk)>0
            tissout(iii,jjj,kkk)=hi(nbou(iii,jjj,kkk));
            else
            tissout(iii,jjj,kkk)=nan;    
            end
        end
    end    
end


A = (reshape(tissout,[length(axt),length(ayt),length(azt)]));

end

