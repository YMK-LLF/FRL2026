function [WEUP] = calculate_W(S,L,one_p,miu)


numerator1 = S*one_p;
denominator1 = one_p'*S*one_p;
wgmvp=numerator1/denominator1;
w1=(1/L)*S*miu;
numerator2=S*one_p*one_p'*S;
denominator2=one_p'*S*one_p;
w2=-(1/L)*(numerator2/denominator2)*miu;


WEUP=wgmvp+w1+w2;