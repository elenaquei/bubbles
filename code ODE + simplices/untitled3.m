% testing the introduction of delays in from_string_to_polynomial_coef

s1 = '-x1 + R_0e^{-px}l2(1-l2)Delay(x2, tau)+ R_0e^{-px}Delay(x3, tau)*(-l5(x1*x1) + (1-2l2)x1) \n';
s2 = '-px3(-x1 + R_0e^{-px}l2(1-l2)Delay(x2, tau)+ R_0e^{-px}Delay(x3, tau)*(-l5(x1*x1) + (1-2l2)x1) ) + \l3x1^0\n';
s3 = '-l5px3(-x1 + R_0e^{-px}l2(1-l2)Delay(x2, tau)+ R_0e^{-px}Delay(x3, tau)*(-l5(x1*x1) + (1-2l2)x1) ) + \l4x1^0';

string_delay = append(append(s1,s2),s3);

from_string_to_polynomial_coef(string_delay)
