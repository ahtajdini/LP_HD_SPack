


module rSeed

lb1, lm1 = (big"1.92672", big"0.886732")
lb2, lm2= (big"8.334426", big"0.99908")
lb3, lm3 = (big"16.588211", big"1.11351256694")
hb1, hm1= (big"-13.3631627",big"0.967038")
hb2, hm2 = ( big"8.3153553",big"1.00001291")
hb3,hm3= (big"31.99434286205", big"1.0329874")


function seed(D)

    if D < 400
        r1 = lb1+ lm1*D/2
        r2 = lb2+ lm2*D/2
        r3 = lb3+ lm3*D/2
        return BigFloat.([r1,r2,r3])

    else
        r1 = hb1+ hm1*D/2
        r2 = hb2+ hm2*D/2
        r3 = hb3+ hm3*D/2
        return BigFloat.([r1,r2,r3])

    end
end
        
end