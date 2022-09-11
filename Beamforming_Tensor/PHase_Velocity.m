function	eV	=	PHase_Velocity(eA,d2)
    [A1,F_r] = size(eA);
    for f = 1:F_r
    eA(:,f) = eA(:,f)/eA(1,f);               
    end
    for inn = 1:F_r
         x_r = angle(eA(:,inn));
        fai_r = zeros(A1,1);
        for in = 1:length(x_r)-1
            if abs(x_r(in+1)-x_r(in))<pi
                fai_r(in+1) = fai_r(in)+x_r(in+1)-x_r(in);
            elseif x_r(in+1)-x_r(in)>pi
                fai_r(in+1) = fai_r(in)+x_r(in+1)-x_r(in)-2*pi;
            else
                fai_r(in+1) = fai_r(in)+x_r(in+1)-x_r(in)+2*pi;
            end
        end
        Aaa = [ones(1,size(eA,1));2*pi*d2]';
        BBB = pinv(Aaa)*fai_r;
        xrr = BBB(2);
        eV(inn) = xrr;
    end 
end