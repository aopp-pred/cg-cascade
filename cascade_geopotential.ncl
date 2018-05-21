; cascade_geopotential.ncl

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


function cascade_geopotential(fout:file,hyam_ec_out:float,hybm_ec_out:float,hyai_ec_out:float,hybi_ec_out:float,ref_P0_out:float)

local fout, hyam_ec_out, hybm_ec_out, hyai_ec_out, hybi_ec_out, ref_P0_out, t_data_out, q_mr_data_out, ps_data_out, z_sfc , ec_pres_hybrid_levs, ec_pres_hybrid_hlevs, epsilon, tv_top, tv_bot , tv, z_levels, linlog, lev_dim,  zh_levels, zhdims


begin

  ;==================================================================================
  ; derived variable: Geopotential height

    print("procedure geopotential height")

    ;=====================
    ; load required variables in from file
    t_data_out              = fout->t
    q_mr_data_out           = fout->q
    ps_data_out             = fout->ps
    z_sfc                   = fout->z_sfc
    ec_pres_hybrid_levs     = fout->pressure_f
    ec_pres_hybrid_hlevs    = fout->pressure_h

    ;====================
    ; calculation

    epsilon = 0.622
    tv_top = t_data_out*(1 + q_mr_data_out/epsilon)
    tv_bot = 1 + q_mr_data_out
    tv = tv_top/tv_bot
    copy_VarCoords(t_data_out,tv)

    ; calculate geopotential: ordered top to bottom (consistent with ECMWF)
    z_levels = cz2ccm(ps_data_out,z_sfc,tv,ref_P0_out,hyam_ec_out(::-1),hybm_ec_out(::-1),hyai_ec_out(::-1),hybi_ec_out(::-1))
    copy_VarCoords(tv,z_levels)

    delete([/epsilon,tv_top,tv_bot/])
    delete([/t_data_out,q_mr_data_out,ps_data_out/])

    ;------
    ; Also calculate on half levels - interpolate
    ; (using same routine gives strange levels near ground)

    linlog = 2   ; 1: linear   2: ln(p) interpolation
    lev_dim = 1
    zh_levels = int2p_n(ec_pres_hybrid_levs,z_levels,ec_pres_hybrid_hlevs,linlog,lev_dim)
    zhdims = dimsizes(zh_levels)
    ; z at surface == z at lowest half level
    zh_levels(:,zhdims(1)-1,:,:) = z_sfc/9.80665
    ; need to approx value zh at level 0 (TOA)
    zh_levels(:,0,:,:) = 2.0*z_levels(:,0,:,:)-zh_levels(:,1,:,:)

    zh_levels!0         = "time"
    zh_levels&time      = z_levels&time
    zh_levels!1         = "nlevp1"
    zh_levels&nlevp1    = hyai_ec_out&nlevp1
    zh_levels!2         = "lat"
    zh_levels&lat       = z_levels&lat
    zh_levels!3         = "lon"
    zh_levels&lon       = z_levels&lon
    zh_levels@long_name = "Height - half level"
    zh_levels@units     = "m"

    copy_VarAtts(zh_levels,z_levels)
    z_levels@long_name = "Height - full level"


    ;-- save to file
    add_to_file(fout,z_levels   , "height_f")
    add_to_file(fout,zh_levels  , "height_h")


    delete([/z_levels,zh_levels,z_sfc/])

    return(tv)

end
