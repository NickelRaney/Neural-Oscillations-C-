for a1 = 1:5
    for a2 = 1:5
        for a3 = 1:5
            for b1 = 1:5
                for b2=1:5
                    for b3=1:5
                        
                        x1 = a1*b1;
                        x2 = a1*b2 + a2*b1;
                        x3 = a1*b3 + b1*a3 + a2*b2;
                        x4 = a2*b3 + a3*b2;
                        x5 = a3*b3;
                        x1 = rem(x1, 5);
                        x2 = rem(x2, 5);
                        x3 = rem(x3, 5);
                        x4 = rem(x4, 5);
                        x5 = rem(x5, 5);
                        if x1==1 && x2==1 &&x3==1 &&x4==0 &&x5==3
                            b1
                            b2
                            b3
                            a1
                            a2
                            a3
                        end
                        
                    end
                end
            end
        end
    end
end