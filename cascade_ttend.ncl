; cascade_ttend.ncl

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


procedure cascade_ttend(fout:file)

local fout, t_data_out, u_data_out, v_data_out, ttend_data_out


begin

  ;==================================================================================
  ; derive advective T tendency

    print("procedure advective T tendency")

    ;=====================
    ; load required variables in from file
    t_data_out           = fout->t
    u_data_out           = fout->u
    v_data_out           = fout->v

    ;=====================
    ; calculation

    ttend_data_out = calc_advection(t_data_out,u_data_out,v_data_out)

    ttend_data_out@long_name = "Advective T tendency" ;
    ttend_data_out@units     = "K/s" ;

    ;-- save to file
    add_to_file(fout,ttend_data_out  , "tadv")

    delete([/t_data_out, u_data_out, v_data_out, ttend_data_out/])

end
