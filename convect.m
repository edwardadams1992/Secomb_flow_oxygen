function [al,cv,segc,isegk,isegkk]=convect(isp,nseg,nnv,nnodbc,bcnod,nodout,oxygen,bcp,solutefac,flowfac,nnodfl,nodrank,nodtyp,nodseg,qq,hd,istart,nspoint,qv,q)


isegkk=zeros(1,nseg);
isegk=0;

for i=1:nnv
    for j=1:nnv
        al(i,j)=0;
    end
end


for j=1:nnodbc
    inod=bcnod(j);
    if nodout(inod)==1
        iseg=nodseg(1,inod);
        if oxygen(isp)==1
            
            segc(iseg)=bloodconc(bcp(j,isp)*solutefac(isp),hd(iseg))*qq(iseg)*flowfac; %%deffo eq 3 with some unit conversion factors
            

        else
            segc(iseg)=bcp(j,isp)*solutefac(isp)*qq(iseg)*flowfac;
        end
        isegkk(iseg)=1;
    end
    
    
end

ineg=0;
ihigh=0;

for in=1:nnodfl
    inod=nodrank(in); % nodrank wrong?
    nodt=nodtyp(inod);
    nout=nodout(inod);
    nin=nodt-nout;
    if nodt>1
        sumin=0;
        hdsumin=0;
        fluxsumin=0;
        for ii=nout+1:nodt
            iseg=nodseg(ii,inod);
            sumin=sumin+qq(iseg)*flowfac;
            hdsumin=hdsumin+qq(iseg)*flowfac*hd(iseg);
            fluxsumin=fluxsumin+segc(iseg);
        end
        if oxygen(isp)==1
            [pb,pp]=blood(fluxsumin/sumin,hdsumin/sumin); 
            
        else
            pb=fluxsumin/sumin;
        end
        sumout=0;
        hdsumout=0;
        
        for ii=1:nout %calc out nodes
            iseg=nodseg(ii,inod);
            isegkk(iseg)=1;
            sumout=sumout+qq(iseg)*flowfac;
            hdsumout=hdsumout+qq(iseg)*flowfac*hd(iseg);
            if oxygen(isp)==1
                segc(iseg)=bloodconc(pb,hd(iseg))*qq(iseg)*flowfac; 
                

            else
                segc(iseg)=pb*qq(iseg)*flowfac;
            end
            if q(iseg)>=0
                i=istart(iseg);
            else
                i=istart(iseg)+nspoint(iseg)-1;
            end
            for jj=nout+1:nodt
                jseg=nodseg(jj,inod);
                if q(jseg)>=0
                    j=istart(jseg)+nspoint(jseg)-1;
                else
                    j=istart(jseg);
                end
                al(i,j)=qq(iseg)*flowfac/sumin;
                if oxygen(isp)==1 && nout>1
                    al(i,j)=al(i,j)*bloodconcp(pb,hd(iseg))/bloodconcp(pb,hdsumin/sumin);     
                end
                for k=1:nnv
                    al(i,k)=al(i,k)+al(i,j)*al(j,k);
                end      
            end      
        end
        
        
    end
  
    for ii=1:nout
        iseg=nodseg(ii,inod);
        for jj=0:nspoint(iseg)-1
            if q(iseg)>=0
                i=istart(iseg)+jj;
            else
                i=istart(iseg)+nspoint(iseg)-jj-1;
                
            end
            segc(iseg)=segc(iseg)-qv(i,isp)/2; %eq3 sort of
            cv(i,isp)=segc(iseg)/qq(iseg)/flowfac;


            
            
            segc(iseg)=segc(iseg)-qv(i,isp)/2;
            
            
        end
        for jj=1:nspoint(iseg)-1
            if q(iseg)>=0
                j=istart(iseg)+jj-1;
                i=istart(iseg)+jj;
            else
                j=istart(iseg)+nspoint(iseg)-jj;
                i=istart(iseg)+nspoint(iseg)-jj-1;
            end
            al(i,j)=1;
            for k=1:nnv
                al(i,k)=al(i,k)+al(i,j)*al(j,k);
            end
        end  
    end
    isegk=isegk+nout;

end



for i=1:nnv
    al(i,i)=0.5;
end



end