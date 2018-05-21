; cascade_geostrophic.ncl

;================================================================================
;   Copyright 2018 Hannah M. Christensen
;
;   Licensed under the Apache License, Version 2.0 (the "License");
;   you may not use this file except in compliance with the License.
;   You may obtain a copy of the License at
;
;       http://www.apache.org/licenses/LICENSE-2.0
;
;   Unless required by applicable law or agreed to in writing, software
;   distributed under the License is distributed on an "AS IS" BASIS,
;   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;   See the License for the specific language governing permissions and
;   limitations under the License.
;================================================================================


procedure cascade_geostrophic(fout:file,tv:float,hyam_ec_out:float,hybm_ec_out:float)



local fout,tv,hyam_ec_out,hybm_ec_out,z_levels,ec_pres_hybrid_levs,grav,R_dry,rpi,rday,rsiyea,rsiday,\
      romega,f_cori_vec,f_cori,geopot,smooth_me,gradphi,gradphi_x,gradphi_y,\
      logpres,gradp,gradp_x,gradp_y,press_gradx_out,press_grady_out,vg_out,ug_out,\
      dsizes,dsizes2,tmparray,lsm_tmp,lsm,lsm_perm,d_size,pgx_out_lsm,pgy_out_lsm,\
      is_cyclic,guess_type,nscan,eps,relc, opt


begin

   ;==================================================================================
   ; param 15 "geostrophic winds" : calculate forcing term from
   ; pressure gradient and gradient in phi
   ;

    print("procedure geostrophic winds")

    ;=====================
    ; load required variables in from file
    z_levels                = fout->height_f       ; dependent on cascade_geopotential.ncl
    ec_pres_hybrid_levs     = fout->pressure_f
    lsm                     = fout->lsm

    ;=====================
    ; do calculation

    ; set up constants - following IFS convention
    grav   = 9.80665
    R_dry  = 287.0
    rpi    = 4*atan(1.0);
    rday   = 86400;
    rsiyea = 365.25*rday*2.0*rpi/6.283076;
    rsiday=rday/(1.0+rday/rsiyea);
    romega=2.0*rpi/rsiday;

    f_cori_vec=2*romega*sin(z_levels&lat*rpi/180);
    f_cori=conform_dims(dimsizes(z_levels),f_cori_vec,2)
    delete([/f_cori_vec/])

    ;----
    ; TERM 1 - partial d phi by dx at constant eta
    geopot = grav*z_levels
    ; smooth it - now do this during hor_derivative
    ; bring back in
 ;hmc;   geopot = smth9(geopot,0.5,0.25,False) ; heavy smoothing
    copy_VarCoords(z_levels,geopot)
    copy_VarAtts(z_levels,geopot)
    smooth_me = False; ;HMC get rid of
    gradphi = hor_derivative(geopot,smooth_me)
    gradphi_x = gradphi[0]
    gradphi_y = gradphi[1]
    delete([/gradphi,geopot/])

    ;----
    ; TERM 2 - (1/rho) * partial d p by dx at constant eta
    logpres = log(ec_pres_hybrid_levs)
    ; smooth it
 ;hmc;   logpres = smth9(logpres,0.5,0.25,False) ; heavy smoothing
    copy_VarCoords(ec_pres_hybrid_levs,logpres)
    copy_VarAtts(ec_pres_hybrid_levs,logpres)
    smooth_me = False
    gradp = hor_derivative(logpres,smooth_me)
    gradp_x = gradp[0]
    gradp_y = gradp[1]
    delete([/gradp,logpres/])

    gradp_x = R_dry*gradp_x*tv
    gradp_y = R_dry*gradp_y*tv

    ;----
    ; COMBINE

    press_gradx_out = gradphi_x + gradp_x
    press_grady_out = gradphi_y + gradp_y

    delete([/grav,R_dry,gradphi_x,gradphi_y,gradp_x,gradp_y/])
    delete([/rpi,rday,rsiyea,rsiday,romega,z_levels/])

;    ; HMC HMC HMC
;    ; test whether this is still necessary
;
;    ;=====================================
;    ; over land mask out geostrophic winds and
;    ; interpolate pressure gradients used to calc ug_out and vg_out
;
;    dsizes = dimsizes(lsm)
;    dsizes2 = (/dsizes(0),dsizes(1),dsizes(2),dimsizes(hyam_ec_out)/)
;    tmparray = new(dsizes2,typeof(lsm))
;    delete([/dsizes,dsizes2/])
;
;    lsm_tmp   = conform(tmparray,lsm,(/0,1,2/))
;    lsm_tmp!0 = "time"
;    lsm_tmp!1 = "lat"
;    lsm_tmp!2 = "lon"
;    lsm_tmp!3 = "nlev"
;    lsm_perm  = lsm_tmp(time | :, nlev | :, lat | :, lon | :)
;    delete([/lsm_tmp/])
;
;    lsm_perm = smth9(lsm_perm,0.5,0.25,False) ; heavy smoothing
;    lsm_perm = smth9(lsm_perm,0.5,0.25,False) ; heavy smoothing
;    lsm_perm = smth9(lsm_perm,0.5,0.25,False) ; heavy smoothing
;    d_size = dimsizes(lsm_perm)
;    lsm_perm(:,:,:,0) = 1            ; mask edge:
;    lsm_perm(:,:,0,:) = 1            ; set values to missing
;    lsm_perm(:,:,:,d_size(3)-1) = 1  ;
;    lsm_perm(:,:,d_size(2)-1,:) = 1  ;
;    pgx_out_lsm = where(lsm_perm.ge.0.0003,press_gradx_out@_FillValue,press_gradx_out)
;    pgy_out_lsm = where(lsm_perm.ge.0.0003,press_grady_out@_FillValue,press_grady_out)
;    copy_VarCoords(press_gradx_out,pgx_out_lsm)
;    copy_VarCoords(press_grady_out,pgy_out_lsm)
;    copy_VarAtts(press_gradx_out,pgx_out_lsm)
;    copy_VarAtts(press_grady_out,pgy_out_lsm)
;    delete([/lsm_perm,press_gradx_out,press_grady_out/])
;
;
;    ;======================
;    ; use poisson grid fill to smoothly fill in missing values
;
;    is_cyclic  = False     ; not cyclic data
;    guess_type = 1         ; start with zonal means
;    nscan      = 200       ; no. iterations
;    eps        = 1.e-6     ; tolerance
;    relc       = 0.6       ; relaxation const
;    opt        = 0         ; dummy
;
;    poisson_grid_fill(pgx_out_lsm,is_cyclic,guess_type,nscan,eps,relc,opt)
;    poisson_grid_fill(pgy_out_lsm,is_cyclic,guess_type,nscan,eps,relc,opt)
;    vg_out =  pgx_out_lsm/f_cori
;    ug_out = -pgy_out_lsm/f_cori
;; new
;
    vg_out =  press_gradx_out/f_cori
    ug_out = -press_grady_out/f_cori
;
;; back to normal:
;
    copy_VarCoords(ec_pres_hybrid_levs,ug_out)
    copy_VarCoords(ec_pres_hybrid_levs,vg_out)
    copy_VarAtts(ec_pres_hybrid_levs,ug_out)
    copy_VarAtts(ec_pres_hybrid_levs,vg_out)
    ug_out@long_name = "Geostrophic U wind"
    ug_out@units     = "m/s"
    vg_out@long_name = "Geostrophic V wind"
    vg_out@units     = "m/s"

    add_to_file(fout,ug_out   , "ug")
    add_to_file(fout,vg_out   , "vg")


;    delete([/pgx_out_lsm,pgy_out_lsm/])
    delete([/ug_out,vg_out,f_cori/])

end
