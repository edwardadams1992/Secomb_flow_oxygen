
function p=eval123(slsegdiv,req,x,nnt,nnv,nsp,mainseg,tisspoints,permsolute,diffsolute,fac,axt,ayt,azt,ds,diff1,g0,qt,scos,qv,ax)
req2=req^2;


for isp=1:nsp
    p(isp)=g0(isp);
end

for itp=1:nnt
    
    dist2=((x(1)-axt(tisspoints(1,itp)))^2)+((x(2)-ayt(tisspoints(2,itp)))^2)+((x(3)-azt(tisspoints(3,itp)))^2);
    if dist2<=req2
        gtt=fac*(1.5-0.5*dist2/req2)/req;
    else
        gtt=fac/(sqrt(dist2));
    end
    
    for isp=1:nsp
        if diffsolute(isp)==1
            p(isp)=p(isp)+gtt/diff1(isp)*qt(itp,isp);
        end
    end
    
end

for i=1:nnv
    iseg=mainseg(i);
    for k=1:slsegdiv
       for j=1:3
           y(j)=ax(j,i)+scos(j,iseg)*ds(iseg)*(-0.5+(k-0.5)/slsegdiv);
       end
       dist2=((x(1)-y(1))^2)+((x(2)-y(2))^2)+((x(3)-y(3))^2);
       if dist2<=req2
           gtv=fac*(1.5-0.5*dist2/req2)/req;
       else
           gtv=fac/sqrt(dist2);
       end
       for isp=1:nsp
           if permsolute(isp)==1
               p(isp)=p(isp)+gtv/diff1(isp)*qv(i,isp)/slsegdiv;
           end
       end
       
    end
end

