function intpk4,nsteps=nsteps1,rangemul=rangemul,nmc=nmc,x1=x1,x2=x2,x3=x3

	if n_elements(rangemul) eq 0 then rangemul=1.0
	

	if n_elements(nsteps1) eq 0 then nsteps=200 else nsteps=nsteps1

	if n_elements(nmc) eq 0 then nmc=floor(nsteps/3.)


	nsteps=long(nsteps)
	
	if n_elements(x1) eq 0 then x1=(20.20)/2.
	if n_elements(x2) eq 0 then x2=(80.78)/2.
	if n_elements(x3) eq 0 then x3=(1301.4)/2.

;	print,x1,x2,x3

	k1max=(0.4+16/x1)*1.25*rangemul
        k2max=(0.4+16/x2)*1.25*rangemul
	k3max=(0.4+16/x3)*1.25*rangemul

; sufficient to encompass entire region where sin^2=1/2 not adequate:
	k1maxb=4.*!pi/x1
	k2maxb=4.*!pi/x2
	k3maxb=4.*!pi/x3

	nstep2=nmc
;floor(nsteps/5)

	dk1=(k1max-k1maxb)/(nstep2-1)
	dk2=(k2max-k2maxb)/(nstep2-1)
	dk3=(k3max-k3maxb)/(nstep2-1)

	dk1b=k1maxb/(nsteps-1)
	dk2b=k2maxb/(nsteps-1)
	dk3b=k3maxb/(nsteps-1)

	k1arr=findgen(nsteps+1)
	k2arr=k1arr*dk2b
	k3arr=k1arr*dk3b
	k1arr=k1arr*dk1b

	w1=sin(k1arr*x1)/k1arr/x1
	w1(where(finite(w1) eq 0))=1.

	w2=sin(k2arr*x2)/k2arr/x2
	w2(where(finite(w2) eq 0))=1.

	w3=sin(k3arr*x3)/k3arr/x3
	w3(where(finite(w3) eq 0))=1.

	w1=w1^2
	w2=w2^2
	w3=w3^2


	k1b=exp(alog(k1max/k1maxb)/nstep2*(findgen(nstep2)+1.))*k1maxb
	k1arr=[k1arr,k1b]
	w1=[w1,(0.5)/x1^2/k1b^2]
	k2b=exp(alog(k2max/k2maxb)/nstep2*(findgen(nstep2)+1.))*k2maxb
	k2arr=[k2arr,k2b]
	w2=[w2,(0.5)/x2^2/k2b^2]
	k3b=exp(alog(k3max/k3maxb)/nstep2*(findgen(nstep2)+1.))*k3maxb
	k3arr=[k3arr,k3b]
	w3=[w1,(0.5)/x3^2/k3b^2]

	integral=0.

	nsteps=nsteps+nstep2
;plot,k3arr

	for i=0,nsteps  do begin
		k1=k1arr(i)
		windowf=w1(i)
		dk1=(k1arr((i+1)<nsteps)-k1arr((i-1)>0))/2.
;		if i eq 0 or i eq nsteps then dk1=dk1*2.
		ktot=k1^2

		ind2=((lindgen(nsteps+1,nsteps+1)) mod (nsteps+1))
		dk2=(k2arr((ind2+1)<nsteps)-k2arr((ind2-1)>0))/2.
;		wh=where(ind2 eq 0 or ind2 eq nsteps)
;		dk2(wh)=dk2(wh)*2.
		ktot=ktot+k2arr(ind2)^2
		windowf=windowf*w2(ind2)
		ind2=1.
		
		ind3=(lindgen(nsteps+1,nsteps+1) / (nsteps+1))
		dk3=(k3arr((ind3+1)<nsteps)-k3arr((ind3-1)>0))/2.
;		wh=where(ind3 eq 0 or ind3 eq nsteps)
;		dk3(wh)=dk3(wh)*2.
		ktot=ktot+k3arr(ind3)^2
		windowf=windowf*w3(ind3)
		ind3=1.
		ktot=sqrt(ktot)

		pk=pofk(ktot)
		if i eq 0 then pk(0,0)=0.
		integral=integral+total(double(pk*windowf*dk1*dk2*dk3))
	endfor



	return,integral*1/!pi^3
end


