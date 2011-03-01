function get_mjdstr
; jm08aug11nyu
    get_juldate, jd
return, string(long(jd-2400000L),format='(I5)')
end
