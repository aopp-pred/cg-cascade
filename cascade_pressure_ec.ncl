; cascade_pressure_ec.ncl

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


function cascade_pressure_ec(ps_data_out:float,ref_P0_out:float,hyam_ec_out:float,hybm_ec_out:float,hyai_ec_out:float,hybi_ec_out:float)


local    dsizes,dsizes2,tmparray,ps_data_tmp,ps_data_out,ps_data_tmp2,ref_P0_out_area,ec_pres_hybrid_levs,ec_pres_hybrid_hlevs

begin

  ;=================================================================================
  ;  define pressure on ECMWF full and half levels
  ; 
  ;  1. FULL LEVELS

    dsizes = dimsizes(ps_data_out)
    dsizes2 = (/dsizes(0),dsizes(1),dsizes(2),1/)
    tmparray = new(dsizes2,typeof(ps_data_out))

    ps_data_tmp   = conform(tmparray,ps_data_out,(/0,1,2/))
    ps_data_tmp!0 = "time"
    ps_data_tmp!1 = "lat"
    ps_data_tmp!2 = "lon"
    ps_data_tmp!3 = "nlev"
    ps_data_tmp2 = ps_data_tmp(time | :, nlev | :, lat | :, lon | :)

    ref_P0_out_area = conform(ps_data_tmp2,ref_P0_out,1)

    ec_pres_hybrid_levs = new((/dsizes(0),dimsizes(hyam_ec_out),dsizes(1),dsizes(2)/),typeof(ps_data_out))

    do i_lev = 0,dimsizes(hyam_ec_out)-1
      ec_pres_hybrid_levs(:,i_lev,:,:) = hyam_ec_out(i_lev)*ref_P0_out_area(:,0,:,:)  +  hybm_ec_out(i_lev)*ps_data_tmp2(:,0,:,:)
    end do

    ec_pres_hybrid_levs!0         = "time"
    ec_pres_hybrid_levs&time      = ps_data_out&time
    ec_pres_hybrid_levs!1         = "nlev"
    ec_pres_hybrid_levs&nlev      =  hyam_ec_out&nlev
    ec_pres_hybrid_levs!2         = "lat"
    ec_pres_hybrid_levs&lat       = ps_data_out&lat
    ec_pres_hybrid_levs!3         = "lon"
    ec_pres_hybrid_levs&lon       = ps_data_out&lon
    ec_pres_hybrid_levs@long_name = "Pressure - full level"
    ec_pres_hybrid_levs@units     = "Pa"


    delete([/dsizes,dsizes2,tmparray,ps_data_tmp,ps_data_tmp2,ref_P0_out_area/])

  ;==================================================================================
  ; param 2. calculate pressure on ECMWF hybrid *half* levels
  ;          not needed for interpolation - only needed for SCM

    dsizes = dimsizes(ps_data_out)
    dsizes2 = (/dsizes(0),dsizes(1),dsizes(2),1/)
    tmparray = new(dsizes2,typeof(ps_data_out))

    ps_data_tmp   = conform(tmparray,ps_data_out,(/0,1,2/))
    ps_data_tmp!0 = "time"
    ps_data_tmp!1 = "lat"
    ps_data_tmp!2 = "lon"
    ps_data_tmp!3 = "nlevp1"
    ps_data_tmp2 = ps_data_tmp(time | :, nlevp1 | :, lat | :, lon | :)

    ref_P0_out_area = conform(ps_data_tmp2,ref_P0_out,1)

    ec_pres_hybrid_hlevs = new((/dsizes(0),dimsizes(hyai_ec_out),dsizes(1),dsizes(2)/),typeof(ps_data_out))

    do i_lev = 0,dimsizes(hyai_ec_out)-1
      ec_pres_hybrid_hlevs(:,i_lev,:,:) = hyai_ec_out(i_lev)*ref_P0_out_area(:,0,:,:)  +  hybi_ec_out(i_lev)*ps_data_tmp2(:,0,:,:)
    end do

    ec_pres_hybrid_hlevs!0         = "time"
    ec_pres_hybrid_hlevs&time      = ps_data_out&time
    ec_pres_hybrid_hlevs!1         = "nlevp1"
    ec_pres_hybrid_hlevs&nlevp1    =  hyai_ec_out&nlevp1
    ec_pres_hybrid_hlevs!2         = "lat"
    ec_pres_hybrid_hlevs&lat       = ps_data_out&lat
    ec_pres_hybrid_hlevs!3         = "lon"
    ec_pres_hybrid_hlevs&lon       = ps_data_out&lon
    ec_pres_hybrid_hlevs@long_name = "Pressure - half level"
    ec_pres_hybrid_hlevs@units     = "Pa"


    delete([/dsizes,dsizes2,tmparray,ps_data_tmp,ps_data_tmp2,ref_P0_out_area/])

    return([/ ec_pres_hybrid_levs,ec_pres_hybrid_hlevs /])

end
