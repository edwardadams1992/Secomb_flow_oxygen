function [output]=Amulti(x,nnod,nodelambda,knowntyp,nodtyp,precond,kmat,hmat,nodnod)


input=x;
for inod=1:nnod
    output(inod)=input(inod);
    
    if knowntyp(inod)~=0
        if knowntyp(inod)~=3
            output(inod)=output(inod)+input(nodelambda(inod));
            output(nodelambda(inod))=input(inod);
        end
        
        for i=1:nodtyp(inod)
            currentnode=nodnod(i,inod);

            if knowntyp(currentnode)~=0
                output(inod)=output(inod)+input(currentnode)*hmat(i+1,inod)*precond(inod)*precond(currentnode);
                if knowntyp(currentnode)~=3
                    output(inod)=output(inod)+input(nodelambda(currentnode))*kmat(i+1,inod)*precond(inod)*precond(nodelambda(currentnode));
                end
                if knowntyp(inod)~=3
                    output(nodelambda(inod))=output(nodelambda(inod))+input(currentnode)*kmat(i+1,inod)*precond(nodelambda(inod))*precond(currentnode);
                end
            end
            
        end
        
        
    end
 
end

end