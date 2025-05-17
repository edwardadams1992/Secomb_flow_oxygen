

function mesh=netflow(mesh)
%define nodes, segments
cnode=(mesh.POS)';
nnod=size(cnode,2);
segnodname=mesh.LINES';
nseg=size(segnodname,2);
nodname=1:nnod;
segname=1:nseg;
nodsegm=6;
diam=double(mesh.DIAM);
hd=0.4*ones(1,nseg);


%bc nodes of known direction
bcsc=mesh.allbcs;
nnodbc=length(bcsc);
bctyp=0*ones(1,nnodbc);
bcnodname=bcsc;
bcprfl=mesh.allbcspress;
bchd=0.4*ones(size(bcprfl));

bifpar=[0.964,6.98,-13.29];
cpar=[0.8,-0.075,-11,12];
viscpar=[6,-0.085,3.2,-2.44,-0.06,0.645];
nitmax=35000;
tol=1e-10;
omega=1.95;
nitmax1=1000;
qtol=1e-3;
hdtol=1e-2;
optw=1.1;
optlam=0.5;
constvisc=3;
vplas=1.0466;
mcv=55;
consthd=0.4;
varyviscosity=1;
phaseseparation=1;
facfp=pi*1333/128/0.01*60/(1e6);
mcvcorr=(92/mcv)^0.33333;



%% setuparrays1

ista=zeros(1,nseg);
iend=zeros(1,nseg);
nk=zeros(1,nnod);
nodout=zeros(1,nnod);
nodrank=zeros(1,nnod);
nodtyp=zeros(1,nnod);

nodnod=zeros(nodsegm,nnod);
nodseg=zeros(nodsegm,nnod);

cond=zeros(1,nseg);
qq=zeros(1,nseg);
tau=zeros(1,nseg);
segpress=zeros(1,nseg);
hdold=zeros(1,nseg);
lseg=zeros(1,nseg);
qold=zeros(1,nseg);
segvar=zeros(1,nseg);
nodvar=zeros(1,nnod);
length_weight=zeros(1,nnod);
nodpress=zeros(1,nnod);
%% analyzenet

for inod=1:nnod
    nodtyp(inod)=0;
end

for iseg=1:nseg
    %search for nodes corresponding to thiis segment
    inod=segnodname(1,iseg); %if node names are squential then no search needed
    if inod<=nnod
        if nodname(inod)==inod
            ista(iseg)=inod;
        else
                for inod=1:nnod
                    if nodname(inod)==segnodname(1,iseg)
                       ista(iseg)=inod;
                    end
                end
            
        end
    end
        
       inod=segnodname(2,iseg); %if node names are squential then no search needed
    if inod<=nnod
        if nodname(inod)==inod
            iend(iseg)=inod;
        else
                for inod=1:nnod
                    if nodname(inod)==segnodname(2,iseg)
                       iend(iseg)=inod;
                    end
                end
          
        end
    end
    
    inod1=ista(iseg);
    inod2=iend(iseg);
    nodtyp(inod1)=nodtyp(inod1)+1;
    nodtyp(inod2)=nodtyp(inod2)+1;
    nodseg(nodtyp(inod1),inod1)=iseg;
    nodseg(nodtyp(inod2),inod2)=iseg;
    nodnod(nodtyp(inod1),inod1)=inod2;
    nodnod(nodtyp(inod2),inod2)=inod1;
    
end

for inodbc=1:nnodbc
    bcnod(inodbc)=0;
    for inod=1:nnod
        if nodname(inod)==bcnodname(inodbc)
            bcnod(inodbc)=inod;
        end
      
    end
end

for iseg=1:nseg
    lseg(iseg)=0;
    
    for k=1:3
        lseg(iseg)=lseg(iseg)+(cnode(k,ista(iseg))-cnode(k,iend(iseg)))^2;
    end
       
end
lseg=sqrt(lseg);

%% flow

nodpress=50*ones(size(nodpress));
relax=1;

for iseg=1:nseg
    q(iseg)=0;
end

for niter=1:nitmax1

    disp(niter)
    if rem(niter,5)==0
        relax=0.8*relax;
    end
    
    for iseg=1:nseg
        qold(iseg)=q(iseg);
        hdold(iseg)=hd(iseg);
        visc=(viscor(diam(iseg),hd(iseg)));
        cond(iseg)=facfp*(diam(iseg)^4)/lseg(iseg)/visc;
    end

    [nodpress,nodtyp,maxerr]=solve(nodsegm,nnod,nodtyp,nodseg,nnodbc,bcprfl,cond,bctyp,nodpress,omega,tol,nodnod,bcnod,nitmax);
    
    for iseg=1:nseg
        q(iseg)=(nodpress(ista(iseg))-nodpress(iend(iseg)))*cond(iseg);
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
            
            
            
            
%       nnodfl=73;      
                   
            
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
            
            
            
            %% dishem
            %boundary nodes - set them to const hd 0.4
            isegk=0;
            for bnod=1:nnodbc
                inod=bcnod(bnod);
                
                if nodout(inod)==1
                    iseg=nodseg(1,inod);
                    hd(iseg)=bchd(bnod);  
                    isegk=isegk+1;
                end
            end
            
            %interior nodes
            for in=1:nnodfl
                inod=nodrank(in);
                nodt=nodtyp(inod);
                nout=nodout(inod);
                nin=nodt-nout;
                if nodt>=2
                    for i=1:nodt
                        segs(i)=nodseg(i,inod);
                        flow(i)=abs(q(segs(i)));
                        
                    end
                end
                
                %2seg nodes
                if nodt==2
                    if nout==2
                        hd(segs(2))=bchd(1);
                    end
                    hd(segs(1))=hd(segs(2));
                end
                
                %3seg nodes
                
                if nodt==3
                    if nout==1
                        hd(segs(1))=(flow(2)*hd(segs(2))+flow(3)*hd(segs(3)))/flow(1);
                    end
                    if nout==2
                        hdd=(1-hd(segs(3)))/diam(segs(3));
                        diaquot=(diam(segs(1))/diam(segs(2)))^2;
                        a=bifpar(3)*(diaquot-1)/(diaquot+1)*hdd;
                        b=1+bifpar(2)*hdd;
                        x0=bifpar(1)*hdd;
                        qikdash=(flow(1)/flow(3)-x0)/(1-2*x0);
                        
                        if qikdash<=0
                            hd(segs(1))=0;
                            hd(segs(2))=hd(segs(3))*flow(3)/flow(2);
                        elseif qikdash>=1
                            hd(segs(2))=0;
                            hd(segs(1))=hd(segs(3))*flow(3)/flow(1);
                        else    
                            rbcrat=1/(1+exp(-a-b*log(qikdash/(1-qikdash))));
                            hd(segs(1))=rbcrat*hd(segs(3))*flow(3)/flow(1);
                            hd(segs(2))=(1-rbcrat)*hd(segs(3))*flow(3)/flow(2);
                        end      
                    end
                    
                    if nout==3
                        
                        for i=1:3
                            hd(segs(i))=bchd(1);
                        end
                        
                    end
                end
                
                if nodt>3
                    if nout==nodt
                        for i=1:nodt
                            hd(segs(i))=bchd(1);
                        end
                    else
                        flowsum=0;
                        hq=0;
                        for i=nout+1:nodt
                            flowsum=flowsum+abs(q(segs(i)));
                            hq=hq+hd(segs(i))*abs(q(segs(i)));
                        end
                       if nout>=1
                           for i=1:nout
                               hd(segs(i))=hq/flowsum;
                           end
                       end
                    end
                end
                
                if nodt>=2
                    isegk=isegk+nout;
                end
      
            end %end of dishehm
            %%
            
   maxqerr=0;
   maxhderr=0;
   errsegq=0;
   errseghd=0;
   
   for iseg=1:nseg
       qchange=q(iseg)-qold(iseg);
       hdchange=hd(iseg)-hdold(iseg);
       hd(iseg)=hdold(iseg)+relax*hdchange;
       
       if abs(qchange)>=maxqerr
           maxqerr=abs(qchange);
           errsegq=iseg;
       end
       if abs(hdchange)>=maxhderr
           maxhderr=abs(hdchange);
           errseghd=iseg;
       end
   end
       
    if maxqerr<qtol && maxhderr<hdtol
        disp('converged2')
        disp(maxhderr)
        break;
    else
        disp('not converged')
    end
end 



%%


if relax<0.9999
    
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
                        iseg=nodseg(j,inod);
                        if (inod==iend(iseg)) && (nk(ista(iseg))==0 || q(iseg)<=0)
                        else
                                                    if (inod==ista(iseg)) && (nk(iend(iseg))==0 || q(iseg)>=0)
                                                        
                                                    else
                                                                        nnodfl=nnodfl+1;
                                                                        nk(inod)=1;
                                                                        nodrank(nnodfl)=inod;  
                                                                        flag=1;
                               
                                                        break
                                                    end
                            
                            
                            
                        end

                    end
                

                end
                
            end
            end
            
            
            
            %% dishem
            %boundary nodes - set them to const hd 0.4
            isegk=0;
            for bnod=1:nnodbc
                inod=bcnod(bnod);
                
                if nodout(inod)==1
                    iseg=nodseg(1,inod);
                    hd(iseg)=bchd(bnod);  
                    isegk=isegk+1;
                end
            end
            
            %interior nodes
            for in=1:nnodfl
                inod=nodrank(in);
                nodt=nodtyp(inod);
                nout=nodout(inod);
                nin=nodt-nout;
                if nodt>=2
                    for i=1:nodt
                        segs(i)=nodseg(i,inod);
                        flow(i)=abs(q(segs(i)));
                        
                    end
                end
                
                %2seg nodes
                if nodt==2
                    if nout==2
                        hd(segs(2))=bchd(1);
                    end
                    hd(segs(1))=hd(segs(2));
                end
                
                %3seg nodes
                
                if nodt==3
                    if nout==1
                        hd(segs(1))=(flow(2)*hd(segs(2))+flow(3)*hd(segs(3)))/flow(1);
                    end
                    if nout==2
                        hdd=(1-hd(segs(3)))/diam(segs(3));
                        diaquot=(diam(segs(1))/diam(segs(2)))^2;
                        a=bifpar(3)*(diaquot-1)/(diaquot+1)*hdd;
                        b=1+bifpar(2)*hdd;
                        x0=bifpar(1)*hdd;
                        qikdash=(flow(1)/flow(3)-x0)/(1-2*x0);
                        
                        if qikdash<=0
                            hd(segs(1))=0;
                            hd(segs(2))=hd(segs(3))*flow(3)/flow(2);
                        elseif qikdash>=1
                            hd(segs(2))=0;
                            hd(segs(1))=hd(segs(3))*flow(3)/flow(1);
                        else    
                            rbcrat=1/(1+exp(-a-b*log(qikdash/(1-qikdash))));
                            hd(segs(1))=rbcrat*hd(segs(3))*flow(3)/flow(1);
                            hd(segs(2))=(1-rbcrat)*hd(segs(3))*flow(3)/flow(2);
                        end      
                    end
                    
                    if nout==3
                        
                        for i=1:3
                            hd(segs(i))=bchd(1);
                        end
                        
                    end
                end
                
                if nodt>3
                    if nout==nodt
                        for i=1:nodt
                            hd(segs(i))=bchd(1);
                        end
                    else
                        flowsum=0;
                        hq=0;
                        for i=nout+1:nodt
                            flowsum=flowsum+abs(q(segs(i)));
                            hq=hq+hd(segs(i))*abs(q(segs(i)));
                        end
                       if nout>=1
                           for i=1:nout
                               hd(segs(i))=hq/flowsum;
                           end
                       end
                    end
                end
                
                if nodt>=2
                    isegk=isegk+nout;
                end
      
            end %end of dishehm
end   
    



mesh=struct('POS',mesh.POS,'LINES',mesh.LINES,'DIAM',mesh.DIAM,'allbcs',mesh.allbcs,'allbcs_press',mesh.allbcspress,'Q',q,'hd',hd);
end
