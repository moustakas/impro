function pofk,karr,gamma=gamma, n=n

if n_elements(n) eq 0 then n = 1.
if n_elements(gamma) eq 0 then gamma=0.21

	omega0=0.3
	sigma8 = 0.58*omega0^(-0.47+0.16*omega0)
        h=1.


; determine sigma8 from model


; following Peacock, p. 528

	keff=(0.172+0.011*(alog(gamma/0.34))^2)*1.
	q=keff/(h*gamma)
	tk=alog(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)^2+(5.46*q)^3+(6.71*q)^4)^(-0.25)
	sigma8squn=keff^(3+n)*tk^2

	q=karr/(h*gamma)   ; k/h Gamma
	tk=alog(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)^2+(5.46*q)^3+(6.71*q)^4)^(-0.25)
	delsq=karr^n*tk^2       ; took out factor of k^3
	delsq=delsq*(sigma8)^2/sigma8squn        ;normalize to Borgani et al. sigma 8, omegam=0.3,omegal=0.7

	pofk=2*!Pi^2*delsq

return,pofk
end
