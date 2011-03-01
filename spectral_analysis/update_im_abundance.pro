function update_im_abundance, hiinolog, hii
; jm04jul22uofa 
; de-bugging code!

    linefit = make_linefit_structure(hiinolog)
    strong = im_abundance(linefit,snrcut=1.0)
    copy_struct, strong, hii
    
return, hii
end
