function skew = my_skew(p)
        skew = [0   -p(3)  p(2);
               p(3)  0    -p(1);
              -p(2)  p(1)   0 ] ;
end