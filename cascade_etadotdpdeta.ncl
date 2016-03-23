; cascade_etadotdpdeta.ncl

procedure cascade_etadotdpdeta(fout:file,hyam_ec_out:float,hybm_ec_out:float)


local fout,hyam_ec_out,hybm_ec_out,ps_data_out,w_data_out,ec_pres_hybrid_levs,\
      ec_pres_hybrid_hlevs,u_data_out,v_data_out,dpsdt,dsizes,dsizes2,tmparray,\
      tmp,dpsdt_surf,hybm_dpsdt,padv,etadotdpdeta_data_out,etadotdpdeta_h,nlevp1


begin

  ;==================================================================================
  ; derived forcing variables:  etadotdpdeta

    print("procedure etadotdpdeta")

    ;=====================
    ; load required variables in from file
    ps_data_out          = fout->ps
    w_data_out           = fout->omega
    ec_pres_hybrid_levs  = fout->pressure_f
    ec_pres_hybrid_hlevs = fout->pressure_h
    u_data_out           = fout->u
    v_data_out           = fout->v

    ;=====================
    ; calculation

    dpsdt = center_finite_diff_n(ps_data_out,ps_data_out&time,False,0,0)
    copy_VarCoords(ps_data_out,dpsdt)
    copy_VarAtts(ps_data_out,dpsdt)

    dsizes = dimsizes(dpsdt)
    dsizes2 = (/dsizes(0),dsizes(1),dsizes(2),1/)
    tmparray = new(dsizes2,typeof(dpsdt))

    tmp   = conform(tmparray,dpsdt,(/0,1,2/))
    tmp!0 = "time"
    tmp!1 = "lat"
    tmp!2 = "lon"
    tmp!3 = "nlev"
    dpsdt_surf = tmp(time | :, nlev | :, lat | :, lon | :)

    hybm_dpsdt = new(dimsizes(w_data_out),typeof(dpsdt))
    do i_lev = 0,dimsizes(hybm_ec_out)-1
      hybm_dpsdt(:,i_lev,:,:) = hybm_ec_out(i_lev)*dpsdt_surf(:,0,:,:)
    end do
    delete([/dsizes,dsizes2,tmparray,tmp,dpsdt_surf,dpsdt/])

    padv = calc_advection(ec_pres_hybrid_levs,u_data_out,v_data_out)
    delete([/u_data_out,v_data_out/])

    etadotdpdeta_data_out = w_data_out + padv - hybm_dpsdt
    delete([/w_data_out,padv,hybm_dpsdt/])

    ;=====================
    ; vert interp to hybrid half levs

    linlog = -2   ; ln(p) interpolation
    lev_dim = 1
    etadotdpdeta_h = int2p_n(ec_pres_hybrid_levs,etadotdpdeta_data_out,ec_pres_hybrid_hlevs,linlog,lev_dim)
    delete([/etadotdpdeta_data_out,ec_pres_hybrid_levs/])

    nlevp1 = ec_pres_hybrid_hlevs&nlevp1
    etadotdpdeta_h(:,0,:,:) = 0.0*etadotdpdeta_h(:,1,:,:)
    etadotdpdeta_h(:,dimsizes(nlevp1)-1,:,:) = 0.0*etadotdpdeta_h(:,1,:,:)

    etadotdpdeta_h!0         = "time"
    etadotdpdeta_h&time      = ps_data_out&time
    etadotdpdeta_h!1         = "nlevp1"
    etadotdpdeta_h&nlevp1    = ec_pres_hybrid_hlevs&nlevp1
    etadotdpdeta_h!2         = "lat"
    etadotdpdeta_h&lat       = ec_pres_hybrid_hlevs&lat
    etadotdpdeta_h!3         = "lon"
    etadotdpdeta_h&lon       = ec_pres_hybrid_hlevs&lon

    ;-- save to file
    add_to_file(fout,etadotdpdeta_h  , "etadotdpdeta")

    delete([/etadotdpdeta_h,ps_data_out/])

end
