    
fTest = zeros(150,1);
for ii = 1:150
    if ii < 30 %0 to 30
        fTest(ii) = 1;
    elseif ii < 60 %31 to 60
        cc = (ii - 30)/(60 - 30);
        fTest(ii) = ((1 - cc)*1 + cc*2);
    elseif ii < 90 %61 to 90
        fTest(ii) = 2;
    elseif ii < 120 %91 to 120
        cc = (ii - 90)/(120 - 90);
        fTest(ii) = ((1 - cc)*2 + cc*3);
    else %121 to 150
        fTest(ii) = 3;
    end
end

plot(fTest)