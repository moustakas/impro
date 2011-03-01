function valid_num_arr, arr, integer=integer
     narr=n_elements(arr)
     nval=intarr(narr)
     for i=0l,narr-1 do nval[i]=valid_num(arr[i],integer=integer)
     return,nval
 end


