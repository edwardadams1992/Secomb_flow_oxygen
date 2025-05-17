


function [mesh]=flow_estimate_bc(mesh,bcprfl)
cnode=(mesh.POS)';
nnod=size(cnode,2);
segnodname=mesh.LINES';
nseg=size(segnodname,2);
nodname=1:nnod;
segname=1:nseg;





bcsc=(mesh.allbcs); %undo
% bcsc=[newedge,mesh.voidbcs'];



%bc nodes of known direction
nnodbck=length(bcsc);
bctyp=0*ones(1,nnodbck);
bcnodname=bcsc;

%% input


nodsegm=6;
maxl=150; 
lb=100; %outerbound distance


diam=double(mesh.DIAM);

hd=0.4*ones(1,nseg);

q=zeros(1,nseg);


segtyp=5*ones(1,nseg);

actual_direction=ones(1,nseg);


bifpar=[0.964,6.98,-13.29];
cpar=[0.8,-0.075,-11,12];
viscpar=[6,-0.085,3.2,-2.44,-0.06,0.645];
fahrpr=[1.7,-0.35,-0.6,-0.01];


optw=1.1;
optlam=0.5;
constvisc=3;
vplas=1.0466;
consthd=0.4;
mcv=55;
mcvcorr=(92/mcv)^0.33333;


solvetype=2;
presstarget1=23.6647;
sheartarget1=2*presstarget1;
varytargetshear=1;
varyviscosity=0;
ktaustart=0.002;
ktausteps=9;
maxinsideit=5;
eps=1e-6;
nitmax2=500000;
nitmax3=25;
omega1=1.95;
omega2=1;
diamcrit=8;
known_flow_weight=10;
seed=11111;

%% setuparrays1

nodout=zeros(1,nnod);
nodnod=zeros(nodsegm,nnod);
nodseg=zeros(nodsegm,nnod);
qq=zeros(1,nseg);
tau=zeros(1,nseg);
segpress=zeros(1,nseg);
nodpress=zeros(1,nnod);
cond=zeros(1,nseg);
condsum=zeros(1,nnod);
hfactor2sum=zeros(1,nnod);
sheartarget=zeros(1,nseg);
shearfac=zeros(1,nseg);
hfactor1=zeros(1,nseg);
hfactor2=zeros(1,nseg);
length_weight=zeros(1,nnod);

nodelambda=zeros(1,nnod);
flow_direction=zeros(1,nseg);
segvar=zeros(1,nseg);
nodvar=zeros(1,nnod);
nk=zeros(1,nnod);
nodrank=zeros(1,nnod);
hmat=zeros(nodsegm+1,nnod);
kmat=zeros(nodsegm+1,nnod);


%% flowtest - find zero flow segments and treat them as unknown boundary nodes
nitmax=10000;
tol=1e-12;
omega=1.95;
rand('seed', 123456);

ista=segnodname(1,:);
iend=segnodname(2,:);

[nodtyp,~] = groupcounts(([ista,iend])');
nodtyp=nodtyp'; %how many connections does each node have

for inod=1:nnod
    [r1,c1]=find(ista==(inod));
    [r2,c2]=find(iend==(inod));   
    nods1=[iend(c1),ista(c2)];
    nodnod(1:nodtyp(inod),inod)=nods1';         
end

%if node is boundary (1 connected point=nodtyp=1) then give random
%pressure, else give preset value
for inod=1:nnod
    
    if nodtyp(inod)==1
        
        nodpress(inod)=randi(32767)*100/(32767);
        
    else
        nodpress(inod)=50;
    end
    
end


 
for niter=1:nitmax
    maxerr=0;
    for inod=1:nnod
        if nodtyp(inod)>1
            pcondsum=0;
            
            for i=1:nodtyp(inod)
                pcondsum=pcondsum+nodpress(nodnod(i,inod))/nodtyp(inod);
      
            end
            press1=omega*(pcondsum-nodpress(inod));
            nodpress(inod)=nodpress(inod)+press1;
            
            if abs(press1)>=maxerr
                maxerr=abs(press1);
                errnode=inod;
            end
            
        end
   
    end
   if maxerr<tol
       break;
   end
end



nnflow=0;
iseg1=0;

for iseg=1:nseg
    if abs(nodpress(ista(iseg))-nodpress(iend(iseg)))<1e-6
        nnflow=nnflow+1;
        inod1=ista(iseg);
        inod2=iend(iseg);
        nodtyp(inod1)=nodtyp(inod1)-1;
        nodtyp(inod2)=nodtyp(inod2)-1;
        
    else
        iseg1=iseg1+1;
        segname(iseg1)=segname(iseg);
        segtyp(iseg1)=segtyp(iseg);
        segnodname(1,iseg1)=segnodname(1,iseg);
        segnodname(2,iseg1)=segnodname(2,iseg);
        diam(iseg1)=diam(iseg);
        q(iseg1)=q(iseg);
        hd(iseg1)=hd(iseg);
        
    end
    
end


nseg=nseg-nnflow;

inod1=0;
nnod0=0;

for inod=1:nnod
    if nodtyp(inod)==0
        nnod0=nnod0+1; 
        
    else
        inod1=inod1+1;
        nodname(inod1)=nodname(inod);
        
        for i=1:3
            cnode(i,inod1)=cnode(i,inod);
        end
      
    end
   
end

nnod=nnod-nnod0;


%% analyzenet

ista=segnodname(1,:);
iend=segnodname(2,:);

[nodtyp,GR] = groupcounts(([ista,iend])');
nodtyp=nodtyp'; %how many connections does each node have

%nodnod - nodes connected to node(col)
%nodseg - segment names for the asscioated node(col)
for inod=1:nnod
    [r1,c1]=find(ista==(inod));
    [r2,c2]=find(iend==(inod));   
    nods1=[iend(c1),ista(c2)];
    nods2=[(c1),(c2)];
    nodnod(1:nodtyp(inod),inod)=nods1';  
    nodseg(1:nodtyp(inod),inod)=nods2';

       
end

%calc lengths of all segments
lseg=sqrt(((cnode(1,ista)-cnode(1,iend)).^2)+((cnode(2,ista)-cnode(2,iend)).^2)+((cnode(3,ista)-cnode(3,iend)).^2));

totallength=sum(lseg);


%weight factors by length - sum over all connected nodes multupled by 0.5
for inod=1:nnod
length_weight(inod)=sum(0.5*lseg(nodseg(1:nodtyp(inod),inod)));
end
totallength1=sum(length_weight);

%classify
numberknownpress=0;

nnodtyp0=sum(nodtyp==0);
nnodbc=sum(nodtyp==1);
numberunknown=sum(nodtyp==1);
bcnod=find(nodtyp==1);
bchd=consthd*ones(size(bcnod));

knowntyp=ones(1,nnod);
knowntyp(nodtyp==1)=3; %boundary nodes =3,i.e intially unknown

lambda=zeros(1,nnod);
nodeinflow=zeros(1,nnod);

for inod=1:nnod
    
nodpress(inod)=presstarget1+randi(32767)*10/32767-5; 

end

% nodpress=[30.888150,28.351756,32.216010,32.576128,28.484512,26.252693,32.087222,35.306009]; %for test - rand works diff in c++

known_flow_direction=zeros(1,nseg);
hd=consthd*ones(1,nseg);

for inodbc=1:nnodbck
    for inod=1:nnod
        if nodname(inod)==bcnodname(inodbc) %if node is bc node
            if bctyp(inodbc)==0 %pressure bc
                knowntyp(inod)=0;
                nodpress(inod)=bcprfl(inodbc);
                numberknownpress=numberknownpress+1;
                numberunknown=numberunknown-1;
                
            end
            if bctyp(inodbc)==2 %flow bc
                
                knowntyp(inod)=2;
                nodeinflow(inod)=bcprfl(inodbc);
                numberunknown=numberunknown-1;
            end
            
            if bctyp(inodbc)==-2 %known direction
                knowntyp(inod)=3;
                iseg=nodseg(1,inod);
                
                if bcprfl(inodbc)>0 && inod==ista(iseg) %if the bcnode is +ve and the start of the segment=outlfow
                    known_flow_direction(iseg)=1;
                elseif bcprfl(inodbc)<0 && inod==iend(iseg) %if the bcnode is -ve and the end of the segment=outflow
                    known_flow_direction(iseg)=1;
                else
                    known_flow_direction(iseg)=-1; %else must be an inflow node - flow moving from the node
                end
                
                
            end
            
            
        end
    end    
end
 
%knowntyp: 0=known pressure, 1=interior node, 2=known flow, 3=known direction, unknown valur
%known_flow_direction: 1=inflow node, -1=outflow node 
matrixdim=2*nnod-numberunknown-numberknownpress;

counter=nnod;
for inod=1:nnod
    if knowntyp(inod)==1 || knowntyp(inod)==2
        counter=counter+1;
        nodelambda(inod)=counter;
    end
    
end

bvector=zeros(1,matrixdim);
xvector=zeros(1,matrixdim);
precond=zeros(1,matrixdim);

%%
kappa=1;

for i=1:ktausteps+3
    
    if i==1 
      ktau=ktaustart;
    elseif i<=ktausteps
        ktau=ktau*2;
    end
    
    insideit=1;
    disp(ktau)
    while true
        
        if i==1 && insideit==1
        rand('seed', seed);
        
        for iseg=1:nseg
            
            if known_flow_direction(iseg)==0
                
                if randi(2)==1
                    flow_direction(iseg)=1;
                else
                    flow_direction(iseg)=-1;
                end
                
            else
                flow_direction(iseg)=known_flow_direction(iseg);
                
            end
%          flow_direction=[-1,-1,1,-1,1,-1,-1];   
        end
        end
        
        
        [q,kappa,lambda,nodpress,jj]=flow1(nseg,varyviscosity,constvisc,diam,lseg,varytargetshear,nodpress,ista,iend,flow_direction,known_flow_direction,ktau,kappa,solvetype,nnod,length_weight,hmat,kmat,nodtyp,nodseg,knowntyp,precond,nodelambda,presstarget1,nodeinflow,matrixdim,nodnod,lambda,nitmax2,totallength,eps,known_flow_weight);
jjstore(i)=jj;
        numdirectionchange=0;

        for iseg=1:nseg
            qf=q(iseg)*flow_direction(iseg);
            if qf<0
                if known_flow_direction(iseg)==0
                    flow_direction(iseg)=-flow_direction(iseg);
                    numdirectionchange=numdirectionchange+1;
                else
                    disp('uh oh!')
                end
            end
        end
        
        
        
        insideit=insideit+1; 
        

       if insideit>maxinsideit  || numdirectionchange==0
           break;  
       end
 
    end
       
end
% 
%% anaylise results
numberzeroflows=0;
numberreversedflow=0;
for iseg=1:nseg
segpress(iseg)=(nodpress(ista(iseg))+nodpress(iend(iseg)))/2;
qq(iseg)=abs(q(iseg));
    if qq(iseg)<1.e-4
        segtyp(iseg)=0;
        numberzeroflows=numberzeroflows+1;
    end
end

for iseg=1:nseg
    if q(iseg)<0
    flow_direction(iseg)=-1;
    else
    flow_direction(iseg)=1;
    end
    
    if flow_direction(iseg)*actual_direction(iseg)==1
    numberreversedflow=numberreversedflow+1;
    end

end

%%
pcaps=0;
totallength1=0;

for iseg=1:nseg
if qq(iseg)>0

    if diam(iseg)<=diamcrit
        pcaps=pcaps+segpress(iseg)*lseg(iseg);
        totallength1=totallength1+lseg(iseg);   
    end


end
end
pcaps=pcaps/totallength1;

for iseg=1:nseg
segvar(iseg)=0;
end

for inodbc=1:nnodbc
inod=bcnod(inodbc);
iseg=nodseg(1,inod);
if ista(iseg)==inod && q(iseg)>0
outflow=0;
elseif iend(iseg)==inod && q(iseg)<0
outflow=0;
else
outflow=1;
end

if diam(iseg)<diamcrit
    if outflow
        pressclass=6;
    else
        pressclass=7;
    end
else
    if segpress(iseg)>pcaps
        if outflow
            pressclass=8;
            else
            pressclass=9;
        end
    else
        if outflow
            pressclass=4;
        else
            pressclass=5;
        end
    end

end
knowntyp(inod)=pressclass;
end




%%


bcidx=knowntyp>1;

newbcs=find(bcidx);
newbcs_press=nodpress(bcidx);

mesh.allbcs=newbcs;
mesh.allbcspress=newbcs_press;


end
